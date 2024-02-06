#ifndef GLOBAL_SPACE_H
#define GLOBAL_SPACE_H

#include <unordered_map>
#include <map>
#include <vector>
#include <deque>
#include <queue>
#include <forward_list>
#include <algorithm>

#include <iostream>
#include <ostream>
#include <cassert>

#include "config.h"

#ifdef CONFIG_PARALLEL
#include "tbb/concurrent_hash_map.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"
#endif

#include "problem.hpp"
#include "clock.hpp"

#include "global/state.hpp"

namespace NP {

	namespace Global {

		template<class Time> class State_space
		{
			public:

			typedef Scheduling_problem<Time> Problem;
			typedef typename Scheduling_problem<Time>::Workload Workload;
			typedef Schedule_state<Time> State;

			static State_space explore(
					const Problem& prob,
					const Analysis_options& opts)
			{
				// doesn't yet support exploration after deadline miss
				assert(opts.early_exit);

				auto s = State_space(prob.jobs, prob.dag, prob.num_processors, opts.timeout,
				                     opts.max_depth, opts.num_buckets);
				s.be_naive = opts.be_naive;
				s.cpu_time.start();
				s.explore();
				s.cpu_time.stop();
				return s;

			}

			// convenience interface for tests
			static State_space explore_naively(
				const Workload& jobs,
				unsigned int num_cpus)
			{
				Problem p{jobs, num_cpus};
				Analysis_options o;
				o.be_naive = true;
				return explore(p, o);
			}

			// convenience interface for tests
			static State_space explore(
				const Workload& jobs,
				unsigned int num_cpus)
			{
				Problem p{jobs, num_cpus};
				Analysis_options o;
				return explore(p, o);
			}

			Interval<Time> get_finish_times(const Job<Time>& j) const
			{
				auto rbounds = rta.find(j.get_id());
				if (rbounds == rta.end()) {
					return Interval<Time>{0, Time_model::constants<Time>::infinity()};
				} else {
					return rbounds->second;
				}
			}

			bool is_schedulable() const
			{
				return !aborted;
			}

			bool was_timed_out() const
			{
				return timed_out;
			}

			unsigned long number_of_states() const
			{
				return num_states;
			}

			unsigned long number_of_edges() const
			{
				return num_edges;
			}

			unsigned long max_exploration_front_width() const
			{
				return width;
			}

			double get_cpu_time() const
			{
				return cpu_time;
			}

			typedef std::deque<State> States;

#ifdef CONFIG_PARALLEL
			typedef tbb::enumerable_thread_specific< States > Split_states;
			typedef std::deque<Split_states> States_storage;
#else
			typedef std::deque< States > States_storage;
#endif

#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH

			struct Edge {
				const Job<Time>* scheduled;
				const State* source;
				const State* target;
				const Interval<Time> finish_range;
				const Interval<Time> start_range;
				const bool segment;
				const Interval<Time> segment_remaining;

				Edge(const Job<Time>* s, const State* src, const State* tgt, const Interval<Time>& sr,
				     const Interval<Time>& fr, const bool seg = false, const Interval<Time>& seg_rem = Interval<Time>{0, 0})
				: scheduled(s)
				, source(src)
				, target(tgt)
				, finish_range(fr)
				, start_range(sr)
				, segment(seg)
				, segment_remaining(seg_rem)
				{
				}

				bool deadline_miss_possible() const
				{
					return scheduled->exceeds_deadline(finish_range.upto());
				}

				Time earliest_finish_time() const
				{
					return finish_range.from();
				}

				Time latest_finish_time() const
				{
					return finish_range.upto();
				}

				Time earliest_start_time() const
				{
//					return finish_range.from() - scheduled->least_cost();
					return start_range.from();
				}

				Time latest_start_time() const
				{
					return start_range.upto();
				}

				bool is_segment() const
				{
					return segment;
				}

				Interval<Time> get_segment_remaining() const
				{
					return segment_remaining;
				}

			};

			const std::deque<Edge>& get_edges() const
			{
				return edges;
			}

			const States_storage& get_states() const
			{
				return states_storage;
			}

#endif
			private:

			typedef State* State_ref;
			typedef typename std::forward_list<State_ref> State_refs;

#ifdef CONFIG_PARALLEL
			typedef tbb::concurrent_hash_map<std::pair<hash_value_t,hash_value_t>, State_refs> States_map;
			typedef typename States_map::accessor States_map_accessor;
#else
			typedef std::unordered_map<std::pair<hash_value_t,hash_value_t>, State_refs> States_map;
#endif

			typedef const Job<Time>* Job_ref;
			typedef std::multimap<Time, Job_ref> By_time_map;

			typedef std::deque<State_ref> Todo_queue;

			typedef Interval_lookup_table<Time, Job<Time>, Job<Time>::scheduling_window> Jobs_lut;

			typedef std::unordered_map<JobID, Interval<Time> > Response_times;

#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
			std::deque<Edge> edges;
#endif

			Response_times rta;

#ifdef CONFIG_PARALLEL
			tbb::enumerable_thread_specific<Response_times> partial_rta;
#endif

			bool aborted;
			bool timed_out;

			const unsigned int max_depth;

			bool be_naive;

			const Workload& jobs;

			// not touched after initialization
			Jobs_lut _jobs_by_win;
			By_time_map _jobs_by_latest_arrival;
			By_time_map _jobs_by_earliest_arrival;
			By_time_map _jobs_by_deadline;
			std::vector<Job_precedence_set> _predecessors;

			// use these const references to ensure read-only access
			const Jobs_lut& jobs_by_win;
			const By_time_map& jobs_by_latest_arrival;
			const By_time_map& jobs_by_earliest_arrival;
			const By_time_map& jobs_by_deadline;
			const std::vector<Job_precedence_set>& predecessors;

			States_storage states_storage;

			States_map states_by_key;
			// updated only by main thread
			unsigned long num_states, width;
			unsigned long current_job_count;
			unsigned long num_edges;

#ifdef CONFIG_PARALLEL
			tbb::enumerable_thread_specific<unsigned long> edge_counter;
#endif
			Processor_clock cpu_time;
			const double timeout;

			const unsigned int num_cpus;

			State_space(const Workload& jobs,
			            const Precedence_constraints &dag_edges,
			            unsigned int num_cpus,
			            double max_cpu_time = 0,
			            unsigned int max_depth = 0,
			            std::size_t num_buckets = 1000)
			: _jobs_by_win(Interval<Time>{0, max_deadline(jobs)},
			               max_deadline(jobs) / num_buckets)
			, jobs(jobs)
			, aborted(false)
			, timed_out(false)
			, be_naive(false)
			, timeout(max_cpu_time)
			, max_depth(max_depth)
			, num_states(0)
			, num_edges(0)
			, width(0)
			, current_job_count(0)
			, num_cpus(num_cpus)
			, jobs_by_latest_arrival(_jobs_by_latest_arrival)
			, jobs_by_earliest_arrival(_jobs_by_earliest_arrival)
			, jobs_by_deadline(_jobs_by_deadline)
			, jobs_by_win(_jobs_by_win)
			, _predecessors(jobs.size())
			, predecessors(_predecessors)
			{
				for (const Job<Time>& j : jobs) {
					_jobs_by_latest_arrival.insert({j.latest_arrival(), &j});
					_jobs_by_earliest_arrival.insert({j.earliest_arrival(), &j});
					_jobs_by_deadline.insert({j.get_deadline(), &j});
					_jobs_by_win.insert(j);
				}

				for (auto e : dag_edges) {
					const Job<Time>& from = lookup<Time>(jobs, e.first);
					const Job<Time>& to   = lookup<Time>(jobs, e.second);
					_predecessors[index_of(to)].push_back(index_of(from));
				}
			}

			private:

			void count_edge()
			{
#ifdef CONFIG_PARALLEL
				edge_counter.local()++;
#else
				num_edges++;
#endif
			}

			static Time max_deadline(const Workload &jobs)
			{
				Time dl = 0;
				for (const auto& j : jobs)
					dl = std::max(dl, j.get_deadline());
				return dl;
			}

			void update_finish_times(Response_times& r, const JobID& id,
			                         Interval<Time> range)
			{
				auto rbounds = r.find(id);
				if (rbounds == r.end()) {
					r.emplace(id, range);
				} else {
					rbounds->second |= range;
				}
				DM("RTA " << id << ": " << r.find(id)->second << std::endl);
			}

			void update_finish_times(
				Response_times& r, const Job<Time>& j, Interval<Time> range)
			{
				update_finish_times(r, j.get_id(), range);
				if (j.exceeds_deadline(range.upto())) {
					DM("*** Deadline miss: " << j << std::endl);
					aborted = true;
				}
			}

			void update_finish_times(const Job<Time>& j, Interval<Time> range)
			{
				Response_times& r =
#ifdef CONFIG_PARALLEL
					partial_rta.local();
#else
					rta;
#endif
				update_finish_times(r, j, range);
			}


			std::size_t index_of(const Job<Time>& j) const
			{
				return (std::size_t) (&j - &(jobs[0]));
			}

			const Job_precedence_set& predecessors_of(const Job<Time>& j) const
			{
				return predecessors[index_of(j)];
			}

			void check_for_deadline_misses(const State& old_s, const State& new_s)
			{
				auto check_from = old_s.core_availability().min();
				auto earliest   = new_s.core_availability().min();

				// check if we skipped any jobs that are now guaranteed
				// to miss their deadline
				for (auto it = jobs_by_deadline.lower_bound(check_from);
				     it != jobs_by_deadline.end(); it++) {
					const Job<Time>& j = *(it->second);
					if (j.get_deadline() < earliest) {
						if (unfinished(new_s, j)) {
							DM("deadline miss: " << new_s << " -> " << j << std::endl);
							// This job is still incomplete but has no chance
							// of being scheduled before its deadline anymore.
							// Abort.
							aborted = true;
							// create a dummy state for explanation purposes
							auto frange = new_s.core_availability() + j.get_cost();
							auto srange = frange - j.get_cost();
							const State& next =
								new_state(new_s, index_of(j), predecessors_of(j),
								          frange, frange, j.get_key());
							// update response times
							update_finish_times(j, frange);
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
							edges.emplace_back(&j, &new_s, &next, srange, frange);
#endif
							count_edge();
							break;
						}
					} else
						// deadlines now after the next earliest finish time
						break;
				}
			}

			void make_initial_state()
			{
				// construct initial state
				states_storage.emplace_back();
				new_state(num_cpus);
			}

			States& states()
			{
#ifdef CONFIG_PARALLEL
				return states_storage.back().local();
#else
				return states_storage.back();
#endif
			}

			template <typename... Args>
			State_ref alloc_state(Args&&... args)
			{
				states().emplace_back(std::forward<Args>(args)...);
				State_ref s = &(*(--states().end()));;

				// make sure we didn't screw up...
				auto njobs = s->number_of_scheduled_jobs();
//				assert (
//					(!njobs && num_states == 0) // initial state
////				    || (njobs == current_job_count + 1) // normal State
//				    || (njobs == current_job_count + 2 && aborted) // deadline miss
//				);

				return s;
			}

			void dealloc_state(State_ref s)
			{
				assert(&(*(--states().end())) == s);
				states().pop_back();
			}

			template <typename... Args>
			State& new_state(Args&&... args)
			{
				return *alloc_state(std::forward<Args>(args)...);
			}

			template <typename... Args>
			State& new_or_merged_state(Args&&... args)
			{
				State_ref s_ref = alloc_state(std::forward<Args>(args)...);

				// try to merge the new state into an existing state
				State_ref s = merge_or_cache(s_ref);
				if (s != s_ref) {
					// great, we merged!
					// clean up the just-created state that we no longer need
					dealloc_state(s_ref);
				}
				return *s;
			}

#ifdef CONFIG_PARALLEL

			// make state available for fast lookup
			void insert_cache_state(States_map_accessor &acc, State_ref s)
			{
				assert(!acc.empty());

				State_refs& list = acc->second;
				list.push_front(s);
			}

			void insert_cache_state( State_ref s){
				States_map_accessor acc;
				if (states_by_key.find(acc, s->get_complete_key())) {
					insert_cache_state(acc, s);
				}else{
					assert(false);
				}
			}

			// returns true if state was merged
			State_ref merge_or_cache(State_ref s)
			{
				States_map_accessor acc;

				while (true) {
					// check if key exists
					if (states_by_key.find(acc, s->get_complete_key())) {
						for (State_ref other : acc->second) {
							if(other->try_to_dominate(*s))
								return other;
							else if (other->try_to_merge(*s))
								return other;
						}
						// If we reach here, we failed to dominate or merge, so go ahead
						// and insert it.
						insert_cache_state(acc, s);
						return s;
					// otherwise, key doesn't exist yet, let's try to create it
					} else if (states_by_key.insert(acc, s->get_complete_key())) {
						// We created the list, so go ahead and insert our state.
						insert_cache_state(acc, s);
						return s;
					}
					// if we raced with concurrent creation, try again
				}
			}

#else

			void cache_state(State_ref s)
			{
				// create a new list if needed, or lookup if already existing
				auto res = states_by_key.emplace(
					std::make_pair(s->get_complete_key(), State_refs()));

				auto pair_it = res.first;
				State_refs& list = pair_it->second;

				list.push_front(s);
			}


			State_ref merge_or_cache(State_ref s_ref)
			{
				State& s = *s_ref;

				const auto pair_it = states_by_key.find(s.get_complete_key());

				// cannot merge if key doesn't exist
				if (pair_it != states_by_key.end())
					for (State_ref other : pair_it->second) {
						if(other->check_reduction_rule(*s_ref)) {
							if (other->try_to_dominate(*s_ref))
								return other;
							else if (other->try_to_merge(*s_ref))
								return other;
						}
					}
				// if we reach here, we failed to merge
				cache_state(s_ref);
				return s_ref;
			}
#endif

			void check_cpu_timeout()
			{
				if (timeout && get_cpu_time() > timeout) {
					aborted = true;
					timed_out = true;
				}
			}

			void check_depth_abort()
			{
				if (max_depth && current_job_count > max_depth)
					aborted = true;
			}

			bool unfinished(const State& s, const Job<Time>& j) const
			{
				return s.job_incomplete(index_of(j));
			}

			bool ready(const State& s, const Job<Time>& j) const
			{
				return unfinished(s, j) && s.job_ready(predecessors_of(j)) && !s.job_preempted(index_of(j));
			}

			bool all_jobs_scheduled(const State& s) const
			{
				return s.number_of_scheduled_jobs() == jobs.size();
			}

			// assumes j is ready
			Interval<Time> ready_times(const State& s, const Job<Time>& j) const
			{
				Interval<Time> r = j.arrival_window();
				for (auto pred : predecessors_of(j)) {
					Interval<Time> ft{0, 0};
					if (!s.get_finish_times(pred, ft))
						ft = get_finish_times(jobs[pred]);
					r.lower_bound(ft.min());
					r.extend_to(ft.max());
				}
				if(s.job_preempted(index_of(j))) {
					Interval<Time> ft{0, 0};
					ft = s.get_segment_finish_time(index_of(j));
					r.lower_bound(ft.min());
					r.extend_to(ft.max());
				}
				return r;
			}

			// assumes j is ready
			Interval<Time> ready_times(
				const State& s, const Job<Time>& j,
				const Job_precedence_set& disregard) const
			{
				Interval<Time> r = j.arrival_window();
				for (auto pred : predecessors_of(j)) {
					// skip if part of disregard
					if (contains(disregard, pred))
						continue;
					Interval<Time> ft{0, 0};
					if (!s.get_finish_times(pred, ft))
						ft = get_finish_times(jobs[pred]);
					r.lower_bound(ft.min());
					r.extend_to(ft.max());
				}
				if(s.job_preempted(index_of(j))) {
					Interval<Time> ft{0, 0};
					// finish times of the job should be in the state
					ft = s.get_segment_finish_time(index_of(j));
					r.lower_bound(ft.min());
//					r.lower_bound(ft.max());
					r.extend_to(ft.max());
				}
				return r;
			}

			Time latest_ready_time(const State& s, const Job<Time>& j) const
			{
				return ready_times(s, j).max();
			}

			Time earliest_ready_time(const State& s, const Job<Time>& j) const
			{
				return ready_times(s, j).min();
			}

			Time latest_ready_time(
				const State& s, Time earliest_ref_ready,
				const Job<Time>& j_hp, const Job<Time>& j_ref) const
			{
				auto rt = ready_times(s, j_hp, predecessors_of(j_ref));
				return std::max(rt.max(), earliest_ref_ready);
			}

			// Find next time by which any job is certainly released.
			// Note that this time may be in the past.
			Time next_higher_prio_job_ready(
				const State& s,
				const Job<Time> &reference_job,
				const Time t_earliest) const
			{
				auto ready_min = earliest_ready_time(s, reference_job);
				Time when = Time_model::constants<Time>::infinity();

				// check everything that overlaps with t_earliest
				for (const Job<Time>& j : jobs_by_win.lookup(t_earliest))
					if (ready(s, j)
					    && j.higher_priority_than(reference_job)) {
						when = std::min(when,
							latest_ready_time(s, ready_min, j, reference_job));
					}

				// let's look at the higher priority preempted jobs
				auto preempted_jobs = s.get_preempted_jobs();
				for (auto it = preempted_jobs.begin(); it != preempted_jobs.end(); it++) {
					const Job<Time>& j = jobs[std::get<0>(*it)];
					if(j.is(reference_job.get_id()))
						continue;
					if (j.higher_priority_than(reference_job)) {
						Interval<Time> ft{0, 0};
						ft = s.get_segment_finish_time(std::get<0>(*it));
						DM("HP -> " << j << " " << ft << std::endl);
						when = std::min(when, ft.max());
					}
				}

				// No point looking in the future when we've already
				// found one in the present.
				if (when <= t_earliest)
					return when;

				// Ok, let's look also in the future.
				for (auto it = jobs_by_latest_arrival
				               .lower_bound(t_earliest);
				     it != jobs_by_latest_arrival.end(); it++) {
					const Job<Time>& j = *(it->second);

					// check if we can stop looking
					if (when < j.latest_arrival())
						break; // yep, nothing can lower 'when' at this point

					// j is not relevant if it is already scheduled or blocked
					if (ready(s, j)
					    && j.higher_priority_than(reference_job)) {
						// does it beat what we've already seen?
						when = std::min(when,
							latest_ready_time(s, ready_min, j, reference_job));
					}
				}

				return when;
			}

			// Find next time by which any job is certainly released.
			// Note that this time may be in the past.
			Time next_job_ready(const State& s, const Time t_earliest) const
			{
				Time when = Time_model::constants<Time>::infinity();

				// check everything that overlaps with t_earliest
				for (const Job<Time>& j : jobs_by_win.lookup(t_earliest))
					if (ready(s, j))
						when = std::min(when, latest_ready_time(s, j));

				// No point looking in the future when we've already
				// found one in the present.
				if (when <= t_earliest)
					return when;

				// Ok, let's look also in the future.
				for (auto it = jobs_by_latest_arrival
				               .lower_bound(t_earliest);
				     it != jobs_by_latest_arrival.end(); it++) {
					const Job<Time>& j = *(it->second);

					// check if we can stop looking
					if (when < j.latest_arrival())
						break; // yep, nothing can lower 'when' at this point

					// j is not relevant if it is already scheduled or blocked
					if (ready(s, j))
						// does it beat what we've already seen?
						when = std::min(when, latest_ready_time(s, j));
				}

				return when;
			}

            // find the next lower or upper bound that a higher priority job after time est can possibly release
			Time possible_preemption(Time est, Time lft, const State &s, const Job<Time> &j) const {
				Time possible_preemption = Time_model::constants<Time>::infinity();
				// first we check lower bounds
				for (auto it = jobs_by_earliest_arrival.lower_bound(est + 1);
					 it != jobs_by_earliest_arrival.upper_bound(lft); it++) {
					const Job<Time> &j_lp = *(it->second);

					// continue if it is the same job
					if(j_lp.is(j.get_id()))
						continue;

					// continue if it is already scheduled
					if (!s.job_incomplete(index_of(j_lp)))
						continue;

					// continue if it is already preempted
					if (s.job_preempted(index_of(j_lp)))
						continue;

					if (j_lp.higher_priority_than(j)) {
						possible_preemption = std::min(possible_preemption, j_lp.earliest_arrival());
						// since we are looking for the next possible preemption, we can stop here
						break;
					}
				}
				// then we check upper bounds
				for (auto it = jobs_by_latest_arrival.lower_bound(est + 1);
					 it != jobs_by_latest_arrival.upper_bound(lft); it++) {
					const Job<Time> &j_hp = *(it->second);

					// continue if it is the same job
					if(j_hp.is(j.get_id()))
						continue;

					// continue if it is already scheduled
					if (!s.job_incomplete(index_of(j_hp)))
						continue;

					// continue if it is already preempted
					if (s.job_preempted(index_of(j_hp)))
						continue;

					if (j_hp.higher_priority_than(j)) {
						possible_preemption = std::min(possible_preemption, j_hp.latest_arrival());
						// since we are looking for the next possible preemption, we can stop here
						break;
					}
				}

				return possible_preemption;
			}

			// check number of possible higher priority job release in an interval
			unsigned int number_of_higher_priority(Time st, Time en, const State &s, const Job<Time> &j) const {
				unsigned int number_of_higher_priority = 0;
				// keep index of higher priority jobs
				std::vector<Job_index> higher_priority_jobs;

				// first we check lower bounds
				for (auto it = jobs_by_earliest_arrival.lower_bound(st+1);
					 it != jobs_by_earliest_arrival.upper_bound(en); it++) {
					const Job<Time> &j_lp = *(it->second);

					// continue if it is the same job
					if(j_lp.is(j.get_id()))
						continue;

					// continue if it is already scheduled
					if (!s.job_incomplete(index_of(j_lp)))
						continue;

					// continue if it is already preempted
					if (s.job_preempted(index_of(j_lp)))
						continue;

					if (j_lp.higher_priority_than(j)) {
						higher_priority_jobs.push_back(index_of(j_lp));
						number_of_higher_priority++;
					}
				}
				// then we check upper bounds
				for (auto it = jobs_by_latest_arrival.lower_bound(st + 1);
					 it != jobs_by_latest_arrival.upper_bound(en); it++) {
					const Job<Time> &j_hp = *(it->second);

					// continue if it is the same job
					if(j_hp.is(j.get_id()))
						continue;

					// continue if it is already scheduled
					if (!s.job_incomplete(index_of(j_hp)))
						continue;

					// continue if it is already preempted
					if (s.job_preempted(index_of(j_hp)))
						continue;

					if (j_hp.higher_priority_than(j) &&
						std::find(higher_priority_jobs.begin(), higher_priority_jobs.end(), index_of(j_hp)) ==
						higher_priority_jobs.end()) {
						number_of_higher_priority++;
					}
				}

				auto preempted_jobs = s.get_preempted_jobs();
				for (auto it = preempted_jobs.begin(); it != preempted_jobs.end(); it++) {
					const Job<Time>& jp = jobs[std::get<0>(*it)];

					// check this finish time if it is in the interval
					if (std::get<2>(*it).min() >= st && std::get<2>(*it).min() <= en) {
						if (jp.higher_priority_than(j)) {
							number_of_higher_priority++;
						}
					}
				}

				return number_of_higher_priority;
			}


			// assumes j is ready
			// NOTE: we don't use Interval<Time> here because the
			//       Interval c'tor sorts its arguments.
			std::pair<Time, Time> start_times(
				const State& s, const Job<Time>& j, Time t_wc, Time &t_high_time, Time &preempt_time) const
			{
				auto rt = ready_times(s, j);
				auto at = s.core_availability();
				Time est = std::max(rt.min(), at.min());

				DM("rt: " << rt << std::endl
				<< "at: " << at << std::endl);

				Time t_preempt = Time_model::constants<Time>::infinity();


				auto t_high = next_higher_prio_job_ready(s, j, at.min());
				DM("t_high: " << t_high << std::endl);
				Time lst    = std::min(t_wc,
					t_high - Time_model::constants<Time>::epsilon());

				// if there is a chance that the job can be dispatched,
				// we calculate the possible preemption point
				if (est <= lst)
					t_preempt = possible_preemption(est,lst + j.get_cost().max(), s, j);

				DM("t_preempt: " << t_preempt << std::endl);
				// now lets see if the job can be preempted
//				if (t_preempt != Time_model::constants<Time>::infinity()) {
//					do {
//						auto cert_avail = s.number_of_certainly_available_cores(t_preempt);
//						auto poss_high = number_of_higher_priority(est, t_preempt, s, j);
//						if (cert_avail - 1 < poss_high) {
//							// the job can be preempted
//							lst = std::min(lst, t_preempt - Time_model::constants<Time>::epsilon());
//							break;
//						} else {
//							// the job cannot be preempted
//							// we have to check the next possible preemption
//							DM("Job cannot be preempted -> look into the next preemption point" << std::endl);
//							t_preempt = possible_preemption(t_preempt, t_preempt + j.get_cost().max(), s, j);
//							DM("Next t_preempt: " << t_preempt << std::endl);
//						}
//					} while (t_preempt < lst + j.get_cost().max());
//				}

				lst = std::min(lst, t_preempt - Time_model::constants<Time>::epsilon());
				preempt_time = t_preempt;
				t_high_time = t_high;

				DM("est: " << est << std::endl);
				DM("lst: " << lst << std::endl);

				return {est, lst};
			}

			bool dispatch(const State& s, const Job<Time>& j, Time t_wc, Interval<Time> exec_time)
			{
				// check if this job has a feasible start-time interval
				Time t_preempt;
				Time t_high;
				auto _st = start_times(s, j, t_wc, t_high, t_preempt);
				if (_st.first > _st.second)
					return false; // nope
				// check if it can also start execution after the possible preemption point
//				if(t_preempt < _st.second)
//					return false; // nope

				Interval<Time> st{_st};

				// yep, job j is a feasible successor in state s

				// compute range of possible finish times
                Interval<Time> ftimes(0, 0);
				ftimes = st + exec_time;

				DM("Assumed finish time: " << ftimes << std::endl);

                Interval<Time> remaining(0, 0);
                if (t_preempt < ftimes.from()) {
					DM("[1] Dispatching segment: " << j << std::endl);
					remaining = {ftimes.from() - t_preempt, ftimes.upto() - t_preempt};
					ftimes = {t_preempt, t_preempt};
				}else if (ftimes.from() < t_preempt && t_preempt < ftimes.until()) {
					DM("[2] Dispatching segment: " << j << std::endl);
					remaining = {0, ftimes.upto() - t_preempt};
					// since job possibly finished at time ftimes.from(),
					// if it is preempted the finish time of the previous segment is t_preempt
					ftimes = {t_preempt, t_preempt};
					// update the finish time
					update_finish_times(j, {ftimes.from(), t_preempt});
				} else {
					// dispatching the whole job
					DM("[3] Dispatching: " << j << std::endl);
				}

                // if we have a leftover, the job is preempted, and we dispatch the first segment
                if (remaining.until() > 0) {
                    // update finish-time
					update_finish_times(j, ftimes);

                    // expand the graph, merging if possible
                    const State& next = be_naive ?
                                        new_state(s, index_of(j), predecessors_of(j),
                                                  st, ftimes, remaining, j.get_key()) :
                                        new_or_merged_state(s, index_of(j), predecessors_of(j),
                                                            st, ftimes, remaining, j.get_key());
                    // make sure we didn't skip any jobs
//                    check_for_deadline_misses(s, next);
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
                    edges.emplace_back(&j, &s, &next,st , ftimes, true, remaining);
#endif
                } else {
                    // update finish-time estimates
                    update_finish_times(j, ftimes);

                    // expand the graph, merging if possible
                    const State &next = be_naive ?
                                        new_state(s, index_of(j), predecessors_of(j),
                                                  st, ftimes, j.get_key()) :
                                        new_or_merged_state(s, index_of(j), predecessors_of(j),
                                                            st, ftimes, j.get_key());
                    // make sure we didn't skip any jobs
//                    check_for_deadline_misses(s, next);
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
                    edges.emplace_back(&j, &s, &next, st, ftimes, false);
#endif
                }



				count_edge();

				return true;
			}

			void explore(const State& s)
			{
				bool found_one = false;

				DM("----" << std::endl);

				// (0) define the window of interest

				// earliest time a core is possibly available
				auto t_min  = s.core_availability().min();
				// latest time some unfinished job is certainly ready
				auto t_job  = next_job_ready(s, t_min);
				// latest time some core is certainly available
				auto t_core = s.core_availability().max();
				// latest time by which a work-conserving scheduler
				// certainly schedules some job
                Time t_wc;
                if (s.has_preempted_jobs()) {
					t_wc = std::max(t_core, t_job);
					// check the first time that the previous segments of a preempted job is completed
					Time min_finish_time = Time_model::constants<Time>::infinity();
					auto preempted_jobs = s.get_preempted_jobs();
					for (auto it = preempted_jobs.begin(); it != preempted_jobs.end(); it++) {
						const Job<Time>& j = jobs[std::get<0>(*it)];
						Interval<Time> ft{0, 0};
						ft = std::get<2>(*it);
						min_finish_time = std::min(min_finish_time, ft.max());
					}
					min_finish_time = std::max(min_finish_time, t_core);
					t_wc = std::min(t_wc, min_finish_time);
				}
                else
				    t_wc = std::max(t_core, t_job);

				DM("Checking state: " << s << std::endl);
				DM("t_min: " << t_min << std::endl
				<< "t_job: " << t_job << std::endl
				<< "t_core: " << t_core << std::endl
				<< "t_wc: " << t_wc << std::endl);

				DM("==== [1] ====" << std::endl);
				// (1) first check jobs that may be already pending
				for (const Job<Time>& j : jobs_by_win.lookup(t_min))
					if (j.earliest_arrival() <= t_min && ready(s, j) && !s.job_preempted(index_of(j)))
                        found_one |= dispatch(s, j, t_wc, j.get_cost());

				DM("==== [2] ====" << std::endl);
				// (2) check jobs that are released only later in the interval
				for (auto it = jobs_by_earliest_arrival.upper_bound(t_min);
					 it != jobs_by_earliest_arrival.end();
					 it++) {
					const Job<Time>& j = *it->second;
					DM(j << " (" << index_of(j) << ")" << std::endl);
					// stop looking once we've left the window of interest
					if (j.earliest_arrival() > t_wc)
						break;

					// Job could be not ready due to precedence constraints
					if (!ready(s, j))
						continue;

                    // if job is previously preempted, we skip it here and we check it later
                    if (s.job_preempted(index_of(j)))
                        continue;

					// Since this job is released in the future, it better
					// be incomplete...
					assert(unfinished(s, j));

					found_one |= dispatch(s, j, t_wc, j.get_cost());
				}

                DM("==== [3] ====" << std::endl);
                // (3) check jobs that are preempted (these jobs are assumed to be certainly released in the past)
                auto preempted_jobs = s.get_preempted_jobs();
				for (auto it = preempted_jobs.begin(); it != preempted_jobs.end(); it++) {
					const Job<Time>& j = jobs[std::get<0>(*it)];
                    found_one |= dispatch(s,  j, t_wc, std::get<1>(*it));
                }

				// check for a dead end
				if (!found_one && !all_jobs_scheduled(s))
					// out of options and we didn't schedule all jobs
					aborted = true;
			}

			// naive: no state merging
			void explore_naively()
			{
				be_naive = true;
				explore();
			}

			void explore() {
				make_initial_state();

				while (current_job_count < jobs.size()) {
					unsigned long n = 0;
#ifdef CONFIG_PARALLEL
					const auto& new_states_part = states_storage.back();
					unsigned int min_jobs = std::numeric_limits<unsigned int>::max();
					for (const States& new_states : new_states_part) {
						for (unsigned int i = 0; i < new_states.size(); i++) {
							const State &s = new_states[i];
							min_jobs = std::min(min_jobs, s.number_of_scheduled_jobs());
						}
					}

					// update current completely scheduled jobs based on the minimum completed jobs in the states
					current_job_count = min_jobs;
#else
					States &exploration_front = states();
//					n = exploration_front.size();
					// update current completely scheduled jobs based on the minimum completed jobs in the states
					auto min_jobs = std::numeric_limits<unsigned int>::max();
					for (const State &s: exploration_front) {
						min_jobs = std::min(min_jobs, s.number_of_scheduled_jobs());
					}
					current_job_count = min_jobs;
#endif

					if (current_job_count == jobs.size())
						break;

					// allocate states space for next depth
					states_storage.emplace_back();


					check_depth_abort();
					check_cpu_timeout();
					if (aborted)
						break;

#ifdef CONFIG_PARALLEL
					// select states with minimum scheduled jobs and explore them
//					std::vector<unsigned int> other_index;
					for (const States& new_states : new_states_part) {
						for (unsigned int i = 0; i < new_states.size(); i++) {
							const State &s = new_states[i];
							unsigned int njobs = s.number_of_scheduled_jobs();
							if (njobs != current_job_count) {
								// copy to next depth
								states().push_back(std::move(s));
//								cache_state(&(*(--states().end())));
								insert_cache_state(&(*(--states().end())));
							}
						}
					}
						// then, explore the states with minimum scheduled jobs in parallel

					// make a local copy of the number of explored states and initialize it to 0
					tbb::enumerable_thread_specific<unsigned long> partial_n;
					for (auto& ni: partial_n){
						ni = 0;
					}
					parallel_for(new_states_part.range(),
						[&] (typename Split_states::const_range_type& r) {
							for (auto it = r.begin(); it != r.end(); it++) {
								const States& new_states = *it;
								auto s = new_states.size();
								tbb::parallel_for(tbb::blocked_range<size_t>(0, s),
									[&] (const tbb::blocked_range<size_t>& r) {
										for (size_t i = r.begin(); i != r.end(); i++)
											if (new_states[i].number_of_scheduled_jobs() == current_job_count){
												explore(new_states[i]);
												partial_n.local()++;
											}
//											explore(new_states[i]);
								});
							}
						});

					// get the number of states explored
					for (auto ni: partial_n){
						n+=ni;
					}
					// keep track of exploration front width
					width = std::max(width, n);

					num_states += n;

#else

					// select states with minimum scheduled jobs and explore them
					std::vector<unsigned int> other_index;

					// Move the states with scheduled jobs more than the minimum to the next depth
					// and keep the rest state index in other_index for later exploration
					for (unsigned int i = 0; i < exploration_front.size(); i++) {
						State& s = exploration_front[i];
						unsigned int njobs = s.number_of_scheduled_jobs();
						if (njobs != min_jobs) {
							// copy to next depth
							states().push_back(std::move(s));
							cache_state(&(*(--states().end())));
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
							// change the state pointer in the edges
							for (auto& e : edges) {
								if (e.source == &s)
									e.source = &states().back();
								if (e.target == &s)
									e.target = &states().back();
							}
#endif
						}
						else
							other_index.push_back(i);
					}

					// then, explore the states with minimum scheduled jobs
					for (auto i : other_index) {
						State& s = exploration_front[i];
						explore(s);
						check_cpu_timeout();
						if (aborted)
							break;
					}

					n = other_index.size();
					// keep track of exploration front width
					width = std::max(width, n);

					num_states += n;


//					for (const State& s : exploration_front) {
//						explore(s);
//						check_cpu_timeout();
//						if (aborted)
//							break;
//					}
#endif

					// clean up the state cache if necessary
					if (!be_naive)
						states_by_key.clear();

					// print number of states
					DM("Number of states: " << num_states << std::endl);
//					std::cout << "states: " << num_states << ", width: " << other_index.size() << ", time: " << get_cpu_time() << std::endl;

#ifdef CONFIG_PARALLEL
					// propagate any updates to the response-time estimates
					for (auto& r : partial_rta)
						for (const auto& elem : r)
							update_finish_times(rta, elem.first, elem.second);
#endif

#ifndef CONFIG_COLLECT_SCHEDULE_GRAPH
					// If we don't need to collect all states, we can remove
					// all those that we are done with, which saves a lot of
					// memory.
#ifdef CONFIG_PARALLEL
					parallel_for(states_storage.front().range(),
						[] (typename Split_states::range_type& r) {
							for (auto it = r.begin(); it != r.end(); it++)
								it->clear();
						});
#endif
					states_storage.pop_front();
#endif
				}


#ifndef CONFIG_COLLECT_SCHEDULE_GRAPH
				// clean out any remaining states
				while (!states_storage.empty()) {
#ifdef CONFIG_PARALLEL
					parallel_for(states_storage.front().range(),
						[] (typename Split_states::range_type& r) {
							for (auto it = r.begin(); it != r.end(); it++)
								it->clear();
						});
#endif
					states_storage.pop_front();
				}
#endif


#ifdef CONFIG_PARALLEL
				for (auto &c : edge_counter)
					num_edges += c;
#endif
			}


#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
			friend std::ostream& operator<< (std::ostream& out,
			                                 const State_space<Time>& space)
			{
					std::map<const Schedule_state<Time>*, unsigned int> state_id;
					unsigned int i = 0;
					std::ostream temp(nullptr);
					out << "digraph {" << std::endl;
#ifdef CONFIG_PARALLEL
					for (const Split_states& states : space.get_states()) {
						for (const Schedule_state<Time>& s : tbb::flattened2d<Split_states>(states)) {
#else
					for (const auto& front : space.get_states()) {
						for (const Schedule_state<Time>& s : front) {
#endif
							if(s.removed())
								continue;
							state_id[&s] = i++;
							out << "\tS" << state_id[&s]
								<< "[label=\"S" << state_id[&s] << ": ";
							s.print_vertex_label(out, space.jobs);
							out << "\"];" << std::endl;
						}
					}
					for (const auto& e : space.get_edges()) {
						out << "\tS" << state_id[e.source]
						    << " -> "
						    << "S" << state_id[e.target]
						    << "[label=\""
						    << "T" << e.scheduled->get_task_id()
						    << " J" << e.scheduled->get_job_id()
						    << "\\nDL=" << e.scheduled->get_deadline()
						    << "\\nES=" << e.earliest_start_time()
 						    << "\\nLS=" << e.latest_start_time()
						    << "\\nEF=" << e.earliest_finish_time()
						    << "\\nLF=" << e.latest_finish_time();
							if (e.is_segment())
								out << "\\nRE: " << e.get_segment_remaining() ;
						    out << "\"";
						if (e.deadline_miss_possible()) {
							out << ",color=Red,fontcolor=Red";
						}
						if (e.is_segment())
							out << ",style=dashed";
						out << ",fontsize=8" << "]"
						    << ";"
						    << std::endl;
						if (e.deadline_miss_possible()) {
							out << "S" << state_id[e.target]
								<< "[color=Red];"
								<< std::endl;
						}
					}
					out << "}" << std::endl;
				return out;
			}
#endif
		};

	}
}

namespace std
{
	template<class Time> struct hash<NP::Global::Schedule_state<Time>>
    {
		std::size_t operator()(NP::Global::Schedule_state<Time> const& s) const
        {
            return s.get_key();
        }
    };
}


#endif
