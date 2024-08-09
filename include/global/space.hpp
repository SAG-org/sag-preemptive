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

namespace Preemptive {

	namespace Global {

		template<class Time>
		class State_space {
		public:

			typedef Scheduling_problem<Time> Problem;
			typedef typename Scheduling_problem<Time>::Workload Workload;
			typedef Schedule_state<Time> State;

			static State_space explore(
					const Problem &prob,
					const Analysis_options &opts) {
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
					const Workload &jobs,
					unsigned int num_cpus) {
				Problem p{jobs, num_cpus};
				Analysis_options o;
				o.be_naive = true;
				return explore(p, o);
			}

			// convenience interface for tests
			static State_space explore(
					const Workload &jobs,
					unsigned int num_cpus) {
				Problem p{jobs, num_cpus};
				Analysis_options o;
				return explore(p, o);
			}

			Interval<Time> get_finish_times(const Job<Time>& j) const
			{
				return Interval<Time>{rta[index_of(j)]};
			}

			bool is_schedulable() const {
				return !aborted;
			}

			bool was_timed_out() const {
				return timed_out;
			}

			unsigned long number_of_states() const {
				return num_states;
			}

			unsigned long number_of_edges() const {
				return num_edges;
			}

			unsigned long max_exploration_front_width() const {
				return width;
			}

			double get_cpu_time() const {
				return cpu_time;
			}

			typedef std::deque<State> States;

#ifdef CONFIG_PARALLEL
			typedef tbb::enumerable_thread_specific<States> Split_states;
			typedef std::deque<Split_states> States_storage;
#else
			typedef std::deque<States> States_storage;
#endif

#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH

			struct Edge {
				const std::vector<const Job<Time>*> scheduled;
				const State* source;
				const State* target;
				const std::vector<Interval<Time>> finish_range;
				const std::vector<Interval<Time>> start_range;
				const bool segment;
				const std::vector<Interval<Time>> segment_remaining;

				Edge(const std::vector<const Job<Time>*> s, const State* src, const State* tgt, const std::vector<Interval<Time>>& sr,
					 const std::vector<Interval<Time>>& fr, const bool seg = false, const std::vector<Interval<Time>>& seg_rem = {Interval<Time>{0, 0}})
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
					bool deadline_miss = false;
					for(int i = 0; i < scheduled.size(); i++) {
						deadline_miss = deadline_miss || scheduled[i]->exceeds_deadline(finish_range[i].upto());
					}
					return deadline_miss;
				}

				std::string earliest_finish_time() const
				{
					std::string finish_time;
					for(int i = 0; i < scheduled.size(); i++) {
						finish_time += std::to_string(finish_range[i].from());
						if(i != scheduled.size() - 1) {
							finish_time += ", ";
						}
					}
					return finish_time;
				}

				std::string latest_finish_time() const
				{
					std::string finish_time;
					for(int i = 0; i < scheduled.size(); i++) {
						finish_time += std::to_string(finish_range[i].upto());
						if(i != scheduled.size() - 1) {
							finish_time += ", ";
						}
					}
					return finish_time;
				}

				std::string earliest_start_time() const
				{
//					return finish_range.from() - scheduled->least_cost();
					std::string start_time;
					for(int i = 0; i < scheduled.size(); i++) {
						start_time += std::to_string(start_range[i].from());
						if(i != scheduled.size() - 1) {
							start_time += ", ";
						}
					}
					return start_time;
				}

				std::string latest_start_time() const
				{
					std::string start_time;
					for(int i = 0; i < scheduled.size(); i++) {
						start_time += std::to_string(start_range[i].upto());
						if(i != scheduled.size() - 1) {
							start_time += ", ";
						}
					}
					return start_time;
				}

				bool is_segment() const
				{
					return segment;
				}

				std::string get_segment_remaining() const
				{
					std::string seg_time;
					for(int i = 0; i < scheduled.size(); i++) {
						seg_time += "[" + std::to_string(segment_remaining[i].from()) + ", " + std::to_string(segment_remaining[i].upto()) + "]";
						if(i != scheduled.size() - 1) {
							seg_time += ", ";
						}
					}
					return seg_time;
				}

				std::string deadline() const
				{
//					return scheduled.back()->get_deadline();
					std::string deadline;
					for(int i = 0; i < scheduled.size(); i++) {
						deadline += std::to_string(scheduled[i]->get_deadline());
						if(i != scheduled.size() - 1) {
							deadline += ", ";
						}
					}
					return deadline;
				}

				unsigned long task_id() const
				{
					return scheduled.back()->get_task_id();
				}

				unsigned long job_id() const
				{
					return scheduled.back()->get_job_id();
				}

				std::string get_label() const {
					std::string label;
					for(int i = 0; i < scheduled.size(); i++) {
						label += "T" +std::to_string(scheduled[i]->get_task_id()) + " J" + std::to_string(scheduled[i]->get_job_id());
						if(i != scheduled.size() - 1) {
							label += ", ";
						}
					}
					return label;
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

			typedef State *State_ref;
			typedef typename std::forward_list<State_ref> State_refs;

#ifdef CONFIG_PARALLEL
			typedef tbb::concurrent_hash_map<std::pair<hash_value_t, hash_value_t>, State_refs> States_map;
			typedef typename States_map::accessor States_map_accessor;
#else
			typedef std::unordered_map<std::pair<hash_value_t, hash_value_t>, State_refs> States_map;
#endif

			typedef const Job<Time> *Job_ref;
			typedef std::multimap<Time, Job_ref> By_time_map;

			typedef std::deque<State_ref> Todo_queue;

			typedef Interval_lookup_table<Time, Job<Time>, Job<Time>::scheduling_window> Jobs_lut;

			// NOTE: we don't use Interval<Time> here because the Interval sorts its arguments.
			typedef std::vector<std::pair<Time, Time>> Response_times;

			typedef std::tuple<Job_ref, Interval<Time>, Time> Job_start;

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

			const Workload &jobs;

			// not touched after initialization
			Jobs_lut _jobs_by_win;
			By_time_map _jobs_by_latest_arrival;
			By_time_map _jobs_by_earliest_arrival;
			By_time_map _jobs_by_deadline;
			std::vector<Job_precedence_set> _predecessors;

			// use these const references to ensure read-only access
			const Jobs_lut &jobs_by_win;
			const By_time_map &jobs_by_latest_arrival;
			const By_time_map &jobs_by_earliest_arrival;
			const By_time_map &jobs_by_deadline;
			const std::vector<Job_precedence_set> &predecessors;

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

			State_space(const Workload &jobs,
						const Precedence_constraints &dag_edges,
						unsigned int num_cpus,
						double max_cpu_time = 0,
						unsigned int max_depth = 0,
						std::size_t num_buckets = 1000)
					: _jobs_by_win(Interval<Time>{0, max_deadline(jobs)},
								   max_deadline(jobs) / num_buckets), jobs(jobs), aborted(false), timed_out(false),
					  be_naive(false), timeout(max_cpu_time), max_depth(max_depth), num_states(0), num_edges(0),
					  width(0), current_job_count(0), num_cpus(num_cpus),
					  jobs_by_latest_arrival(_jobs_by_latest_arrival),
					  jobs_by_earliest_arrival(_jobs_by_earliest_arrival), jobs_by_deadline(_jobs_by_deadline),
					  jobs_by_win(_jobs_by_win), _predecessors(jobs.size()), predecessors(_predecessors),
					  rta(Response_times(jobs.size(), {Time_model::constants<Time>::infinity(), 0}))
#ifdef CONFIG_PARALLEL
			, partial_rta(Response_times(jobs.size(), {Time_model::constants<Time>::infinity(), 0}))
#endif
			{
				for (const Job<Time> &j: jobs) {
					_jobs_by_latest_arrival.insert({j.latest_arrival(), &j});
					_jobs_by_earliest_arrival.insert({j.earliest_arrival(), &j});
					_jobs_by_deadline.insert({j.get_deadline(), &j});
					_jobs_by_win.insert(j);
				}

				for (auto e: dag_edges) {
					const Job<Time> &from = lookup<Time>(jobs, e.first);
					const Job<Time> &to = lookup<Time>(jobs, e.second);
					_predecessors[index_of(to)].push_back(index_of(from));
				}
			}

		private:

			void count_edge() {
#ifdef CONFIG_PARALLEL
				edge_counter.local()++;
#else
				num_edges++;
#endif
			}

			static Time max_deadline(const Workload &jobs) {
				Time dl = 0;
				for (const auto &j: jobs)
					dl = std::max(dl, j.get_deadline());
				return dl;
			}

			void update_finish_times(Response_times &r, const Job_index index,
									 Interval<Time> range) {
				r[index] = std::pair<Time, Time>{std::min(r[index].first, range.from()),
												 std::max(r[index].second, range.upto())};
				DM("RTA " << index << ": [" << r[index].first << ", " << r[index].second << "]" << std::endl);
			}

			void update_finish_times(Response_times &r, const Job_index index,
									 std::pair<Time, Time> range) {
				r[index] = std::pair<Time, Time>{std::min(r[index].first, range.first),
												 std::max(r[index].second, range.second)};
				DM("RTA " << index << ": [" << r[index].first << ", " << r[index].second << "]" << std::endl);
			}

			void update_finish_times(
					Response_times &r, const Job<Time> &j, Interval<Time> range) {
				update_finish_times(r, index_of(j), range);
				if (j.exceeds_deadline(range.upto())) {
					DM("*** Deadline miss: " << j << std::endl);
					aborted = true;
				}
			}

			void update_finish_times(const Job<Time> &j, Interval<Time> range) {
				Response_times &r =
#ifdef CONFIG_PARALLEL
						partial_rta.local();
#else
						rta;
#endif
				update_finish_times(r, j, range);
			}


			std::size_t index_of(const Job<Time> &j) const {
				// make sure that the job is part of the workload
				// and catch the case where the job is not part of the workload,
				// but the user tries to access it anyway
				auto index = (std::size_t) (&j - &(jobs[0]));
				try {
					jobs.at(index);
				} catch (std::out_of_range &e) {
					std::cerr << "Job " << j << " not found in workload." << std::endl;
					std::abort();
				}
				return index;
			}

			const Job_precedence_set &predecessors_of(const Job<Time> &j) const {
				return predecessors[index_of(j)];
			}

			void check_for_deadline_misses(const State &old_s, const State &new_s) {
				auto check_from = old_s.core_availability().min();
				auto earliest = new_s.core_availability().min();

				// check if we skipped any jobs that are now guaranteed
				// to miss their deadline
				for (auto it = jobs_by_deadline.lower_bound(check_from);
					 it != jobs_by_deadline.end(); it++) {
					const Job<Time> &j = *(it->second);
					if (j.get_deadline() < earliest) {
						if (unfinished(new_s, j)) {
							DM("deadline miss: " << new_s << " -> " << j << std::endl);
							// This job is still incomplete but has no chance
							// of being scheduled before its deadline anymore.
							// Abort.
							aborted = true;
							// create a dummy state for explanation purposes
							auto frange = new_s.core_availability() + j.get_cost();
							//auto srange = frange - j.get_cost();
							const State &next =
									new_state(new_s, index_of(j), predecessors_of(j),
											  frange, frange, j.get_key());
							// update response times
							update_finish_times(j, frange);
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
							std::vector<Job_ref> selected_jobs;
							selected_jobs.push_back(&j);
							std::vector<Interval<Time>> finish_times;
							std::vector<Interval<Time>> start_times;
							finish_times.push_back(frange);
							start_times.push_back(new_s.core_availability());
							edges.emplace_back(selected_jobs, &new_s, &next,start_times, finish_times, false);
#endif
							count_edge();
							break;
						}
					} else
						// deadlines now after the next earliest finish time
						break;
				}
			}

			void make_initial_state() {
				// construct initial state
				states_storage.emplace_back();
				new_state(num_cpus);
			}

			States &states() {
#ifdef CONFIG_PARALLEL
				return states_storage.back().local();
#else
				return states_storage.back();
#endif
			}

			template<typename... Args>
			State_ref alloc_state(Args &&... args) {
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

			void dealloc_state(State_ref s) {
				assert(&(*(--states().end())) == s);
				states().pop_back();
			}

			template<typename... Args>
			State &new_state(Args &&... args) {
				return *alloc_state(std::forward<Args>(args)...);
			}

			template<typename... Args>
			State &new_or_merged_state(Args &&... args) {
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
			void insert_cache_state(States_map_accessor &acc, State_ref s) {
				assert(!acc.empty());

				State_refs &list = acc->second;
				list.push_front(s);
			}

			// returns true if state was merged
			State_ref merge_or_cache(State_ref s) {
				States_map_accessor acc;

				while (true) {
					// check if key exists
					if (states_by_key.find(acc, s->get_complete_key())) {
						for (State_ref other: acc->second) {
							if (other->check_reduction_rule(*s)) {
//								if (other->try_to_dominate(*s))
//									return other;
//								else
								if (other->try_to_merge(*s))
									return other;
							}
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

			void cache_state(State_ref s) {
				// create a new list if needed, or lookup if already existing
				auto res = states_by_key.emplace(
						std::make_pair(s->get_complete_key(), State_refs()));

				auto pair_it = res.first;
				State_refs &list = pair_it->second;

				list.push_front(s);
			}


			State_ref merge_or_cache(State_ref s_ref) {
				State &s = *s_ref;

				const auto pair_it = states_by_key.find(s.get_complete_key());

				// cannot merge if key doesn't exist
				if (pair_it != states_by_key.end())
					for (State_ref other: pair_it->second) {
						if (other->check_reduction_rule(*s_ref)) {
//							if (other->try_to_dominate(*s_ref))
//								return other;
//							else
							if (other->try_to_merge(*s_ref))
								return other;
						}
					}
				// if we reach here, we failed to merge
				cache_state(s_ref);
				return s_ref;
			}

#endif

			void check_cpu_timeout() {
				if (timeout && get_cpu_time() > timeout) {
					aborted = true;
					timed_out = true;
				}
			}

			void check_depth_abort() {
				if (max_depth && current_job_count > max_depth)
					aborted = true;
			}

			bool unfinished(const State &s, const Job<Time> &j) const {
				return s.job_incomplete(index_of(j));
			}

			bool ready(const State &s, const Job<Time> &j) const {
				return unfinished(s, j) && s.job_ready(predecessors_of(j)) && !s.job_preempted(index_of(j));
			}

			bool all_jobs_scheduled(const State &s) const {
				return s.number_of_scheduled_jobs() == jobs.size();
			}

			// assumes j is ready
			Interval<Time> ready_times(const State &s, const Job<Time> &j) const {
				Interval<Time> r = j.arrival_window();
				for (auto pred: predecessors_of(j)) {
					Interval<Time> ft{0, 0};
					if (!s.get_finish_times(pred, ft))
						ft = get_finish_times(jobs[pred]);
					r.lower_bound(ft.min());
					r.extend_to(ft.max());
				}
				if (s.job_preempted(index_of(j))) {
					Interval<Time> ft{0, 0};
					ft = s.get_segment_finish_time(index_of(j));
					r.lower_bound(ft.min());
					r.extend_to(ft.max());
				}
				return r;
			}

			// assumes j is ready
			Interval<Time> ready_times(
					const State &s, const Job<Time> &j,
					const Job_precedence_set &disregard) const {
				Interval<Time> r = j.arrival_window();
				for (auto pred: predecessors_of(j)) {
					// skip if part of disregard
					if (contains(disregard, pred))
						continue;
					Interval<Time> ft{0, 0};
					if (!s.get_finish_times(pred, ft))
						ft = get_finish_times(jobs[pred]);
					r.lower_bound(ft.min());
					r.extend_to(ft.max());
				}
				if (s.job_preempted(index_of(j))) {
					Interval<Time> ft{0, 0};
					// finish times of the job should be in the state
					ft = s.get_segment_finish_time(index_of(j));
					r.lower_bound(ft.min());
//					r.lower_bound(ft.max());
					r.extend_to(ft.max());
				}
				return r;
			}

			Time latest_ready_time(const State &s, const Job<Time> &j) const {
				return ready_times(s, j).max();
			}

			Time earliest_ready_time(const State &s, const Job<Time> &j) const {
				return ready_times(s, j).min();
			}

			Time latest_ready_time(
					const State &s, Time earliest_ref_ready,
					const Job<Time> &j_hp, const Job<Time> &j_ref) const {
				auto rt = ready_times(s, j_hp, predecessors_of(j_ref));
				return std::max(rt.max(), earliest_ref_ready);
			}

			// Find the next time by which any job is certainly released.
			// Note that this time may be in the past.
			Time next_higher_prio_job_ready(
					const State &s,
					const Job<Time> &reference_job,
					const Time t_earliest) const {
				auto ready_min = earliest_ready_time(s, reference_job);
				Time when = Time_model::constants<Time>::infinity();

				// check everything that overlaps with t_earliest
				for (const Job<Time> &j: jobs_by_win.lookup(t_earliest))
					if (ready(s, j)
						&& j.higher_priority_than(reference_job)) {
						when = std::min(when,
										latest_ready_time(s, ready_min, j, reference_job));
					}

				// let's look at the higher priority preempted jobs
				auto preempted_jobs = s.get_preempted_jobs();
				for (auto it = preempted_jobs.begin(); it != preempted_jobs.end(); it++) {
					const Job<Time> &j = jobs[std::get<0>(*it)];
					if (j.is(reference_job.get_id()))
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
					const Job<Time> &j = *(it->second);

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

			// Find the next time by which any job is certainly released.
			// Note that this time may be in the past.
			Time next_job_ready(const State &s, const Time t_earliest) const {
				Time when = Time_model::constants<Time>::infinity();

				// check everything that overlaps with t_earliest
				for (const Job<Time> &j: jobs_by_win.lookup(t_earliest))
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
					const Job<Time> &j = *(it->second);

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
				int k = 0;

				// first, let's see how many higher priority jobs already exist in EST
				for (const Job<Time>& j_t : jobs_by_win.lookup(est + Time_model::constants<Time>::epsilon())) {
					// continue if it is already scheduled
					if (!s.job_incomplete(index_of(j_t)))
						continue;

//					if( est < j_t.earliest_arrival())
//						continue;

					// continue if it is already preempted
					if (s.job_preempted(index_of(j_t))) {
						// we have to check the finish time of its previous segment
						// to see if it can block the current job
						Interval<Time> ft = s.get_segment_finish_time(index_of(j_t));
						if (ft.max() < est)
							continue;
					}

					if (j_t.higher_priority_than(j)) {
						k++;
					}

//					if (k >= num_cpus || est < s.core_availability(k).max()) {
//						return est;
//					}
				}

				// then we check lower bounds
				for (auto it = jobs_by_earliest_arrival.lower_bound(est + Time_model::constants<Time>::epsilon());
					 it != jobs_by_earliest_arrival.upper_bound(lft - Time_model::constants<Time>::epsilon()); it++) {
					const Job<Time> &j_lp = *(it->second);

					// continue if it is the same job // GN: that should already be taken care of by th test higher_priority_than(.)
					//if (j_lp.is(j.get_id()))
					//	continue;

					// continue if it is already scheduled
					if (!s.job_incomplete(index_of(j_lp)))
						continue;

					// continue if it is already preempted
					if (s.job_preempted(index_of(j_lp)))
						continue;

					if (j_lp.higher_priority_than(j)) {
						k++;
						if (k >= num_cpus || it->first < s.core_availability(k).max())
						{
							possible_preemption = j_lp.earliest_arrival();
							break;
						}
					}
				}
				// then we check upper bounds
				k = 0;
				for (auto it = jobs_by_latest_arrival.lower_bound(est + Time_model::constants<Time>::epsilon());
					 it != jobs_by_latest_arrival.upper_bound(lft - Time_model::constants<Time>::epsilon()); it++) {
					const Job<Time> &j_hp = *(it->second);

					// continue if it is the same job
					//if (j_hp.is(j.get_id()))
					//	continue;

					// continue if it is already scheduled
					if (!s.job_incomplete(index_of(j_hp)))
						continue;

					// continue if it is already preempted
					if (s.job_preempted(index_of(j_hp)))
						continue;

					if (j_hp.higher_priority_than(j)) {
						k++;
						if (k >= num_cpus || it->first < s.core_availability(k).max())
						{
							possible_preemption = std::min(possible_preemption, j_hp.latest_arrival());
							break;
						}
					}
				}

				return possible_preemption;
			}

			Time calculate_twc(const State& s, const Time t_earliest, const int batch_size) const
			{
				// keeping track when the next 'batch_size' jobs will certainly ready in a non-decreasing ready time order
				std::vector<Time> when(batch_size, Time_model::constants<Time>::infinity());

				// check everything that was certainly released before t_earliest
				Time t_min_rel = t_earliest;
				for (const Job<Time>& j : jobs_by_win.lookup(t_earliest))
				{
					// j is not relevant if it is already scheduled or blocked
					if (ready(s, j)) {
						t_min_rel = std::min(t_min_rel, j.latest_arrival());
					}
				}

				// look in the future only when we've not already
				// found 'batch_size' ready jobs in the present.
//				if (when[batch_size - 1] > t_earliest)
//				{
				int i = 0;
				auto it = jobs_by_latest_arrival.lower_bound(t_min_rel);
				//for (int k = 0; k < batch_size; k++) {
					for (; it != jobs_by_latest_arrival.end(); it++) {
						const Job<Time>& j = *(it->second);
						Time a_max = it->first;

						// check if we can stop looking
						if (when[batch_size - 1] < a_max)
							break; // yep, nothing can lower when 'batch_size' jobs will be ready at this point

						// j is not relevant if it is already scheduled or blocked
						if (ready(s, j)) {
							/*Time rt = latest_ready_time(s, j); //BUG: does not work with DAGs
							if (rt < when[batch_size - 1])
							{
								when[batch_size - 1] = rt;
								std::sort(when.begin(), when.end());
							}*/
							when[i] = a_max;
							i++;
							if (i >= batch_size )
								break;
						}
					}
				//}
				//}

				//if (s.has_preempted_jobs()) {
					// check when the previous segments of each preempted job is finished
					auto preempted_jobs = s.get_preempted_jobs();
					for (auto p = preempted_jobs.begin(); p != preempted_jobs.end(); p++) {
						Interval<Time> ft = std::get<2>(*p);
						if (ft.max() < when[batch_size - 1])
						{
							when[batch_size - 1] = ft.max();
							std::sort(when.begin(), when.end());
						}
					}
				//}

				for (int k = 0; k < batch_size - 1; k++) {
					auto delta_k = std::max(s.core_availability(k).max(), when[k]);
					if (delta_k < s.core_availability(k + 1).min()) {
						return delta_k;
					}
				}
				return std::max(s.core_availability(batch_size - 1).max(), when[batch_size - 1]);
			}

			Time calculate_t_high(
				const State& s,
				const Job<Time>& reference_job,
				const Time t_earliest,
				const Time upper_bound,
				const int batch_size) const
			{
				auto ready_min = earliest_ready_time(s, reference_job);

				// keeping track when the next 'batch_size' higher-priority jobs will certainly ready in a non-decreasing ready time order
				std::vector<Time> when(batch_size, upper_bound);

				// check everything that was certainly released before t_earliest
				Time t_min_rel = t_earliest;
				for (const Job<Time>& j : jobs_by_win.lookup(t_earliest))
				{
					// j is not relevant if it is already scheduled or blocked
					if (ready(s, j)) {
						t_min_rel = std::min(t_min_rel, j.latest_arrival());
					}
				}

				// look in the future only when we've not already
				// found 'batch_size' ready jobs in the present.
//				if (when[batch_size - 1] > t_earliest)
//				{
				int i = 0;
				for (auto it = jobs_by_latest_arrival.lower_bound(t_min_rel); it != jobs_by_latest_arrival.end(); it++) {
					const Job<Time>& j = *(it->second);
					Time a_max = it->first; // latest arrival time

					// check if we can stop looking
					if (when[batch_size - 1] < a_max)
						break; // yep, nothing can lower when 'batch_size' jobs will be ready at this point

					// if we passed the upper bound on the time window we must check, we return the upper bound
					if (a_max >= upper_bound)
						break;

					// j is not relevant if it is already scheduled or blocked
					if (ready(s, j) && j.higher_priority_than(reference_job)) {
						/*Time rt = latest_ready_time(s, ready_min, j, reference_job); //BUG: does not work with DAGs
						if (rt < when[batch_size - 1])
						{
							when[batch_size - 1] = rt;
							std::sort(when.begin(), when.end());
						}*/
						when[i] = a_max;
						i++;
						if (i >= batch_size)
							break;
					}
				}
				//}

				//if (s.has_preempted_jobs()) {
					// check when the previous segments of each preempted job is finished
					auto preempted_jobs = s.get_preempted_jobs();
					for (auto p = preempted_jobs.begin(); p != preempted_jobs.end(); p++) {
						const Job<Time>& j = jobs[std::get<0>(*p)];
						Interval<Time> ft = std::get<2>(*p);
						if (j.higher_priority_than(reference_job) && ft.max() < when[batch_size - 1])
						{
							when[batch_size - 1] = ft.max();
							std::sort(when.begin(), when.end());
						}
					}
				//}

				for (int k = 0; k < batch_size - 1; k++) {
					if (when[k] < s.core_availability(k + 1).min())
						return when[k];
				}
				return when[batch_size - 1];
			}

			unsigned int max_batch_size(const State& s, const Time t_earliest) {

				Time t_min_rel = t_earliest;
				for (const Job<Time>& j : jobs_by_win.lookup(t_earliest))
				{
					// j is not relevant if it is already scheduled or blocked
					if (ready(s, j)) {
						t_min_rel = std::min(t_min_rel, j.earliest_arrival());
					}
				}

				std::vector<Time> when(num_cpus, Time_model::constants<Time>::infinity());
				auto it = jobs_by_earliest_arrival.lower_bound(t_min_rel);
				for (; it != jobs_by_earliest_arrival.end(); it++) {
					const Job<Time>& j = *(it->second);
					// j is not relevant if it is already scheduled or blocked
					if (ready(s, j)) {
						when[0] = it->first;
						break;
					}
				}
				//if (it == jobs_by_earliest_arrival.end())
				//	return 1;

				for (int k = 1; k < num_cpus && it != jobs_by_earliest_arrival.end(); k++) {
					for (it++; it != jobs_by_earliest_arrival.end(); it++) {
						const Job<Time>& j = *(it->second);
						// j is not relevant if it is already scheduled or blocked
						if (ready(s, j)) {
							when[k] = it->first;
							break;
						}
					}
				}

				// check when the previous segments of each preempted job is finished
				auto preempted_jobs = s.get_preempted_jobs();
				for (auto p = preempted_jobs.begin(); p != preempted_jobs.end(); p++) {
					const Job<Time>& j = jobs[std::get<0>(*p)];
					Interval<Time> ft = std::get<2>(*p);
					if (ft.min() < when[num_cpus - 1])
					{
						when[num_cpus - 1] = ft.min();
						std::sort(when.begin(), when.end());
					}
				}

				for (int k = 1; k < num_cpus; k++) {
					if (when[k] < s.core_availability(k).max())
						return k;
				}
				return num_cpus;
			}


			// Helper function to generate combinations
			void combine_helper(const std::vector<Job_start>& set, int start, int k, std::vector<Job_start>& current, std::vector<std::vector<Job_start>>& result, Time threshold) {
				// If the combination is done
				if (k == 0) {
					// if the earliest start time of the last job added is later than the latest start time of any job
					// that is not in the set of selected jobs, than this combination is not legal and we do not add it to the resulting combinations
					assert(start > 0);
					Interval<Time> st = std::get<1>(set[start - 1]);
					for (int i = start; i < set.size(); ++i) {
						Interval<Time> st_other = std::get<1>(set[i]);
						if (st.from() > st_other.max())
							return;
					}
					result.push_back(current);
					return;
				}
				// Try every element to fill the current position
				for (int i = start; i <= set.size() - k; ++i) {
					// if the earliest start time of the next job is later than the latest start time of any job
					// that is not in the set of selected job, than this combination is not legal and we can stop
					Interval<Time> st = std::get<1>(set[i]);
					if (st.from() > threshold)
						break;
					current.push_back(set[i]);
					combine_helper(set, i + 1, k - 1, current, result, threshold);
					current.pop_back();  // backtrack
					// we update the threshold
					threshold = std::min(threshold, st.max());
				}
			}

			// Main function to return combinations of k elements in a vector
			std::vector<std::vector<Job_start>> combine(std::vector<Job_start> &selected_jobs, int k) {
				std::sort(selected_jobs.begin(), selected_jobs.end(), [](const Job_start& a, const Job_start& b) {
					return std::get<1>(a).from() < std::get<1>(b).from();
				});

				// find the maximum size of the combinations that can be generated C(jobs.size(), k)
				int max_size = 1;
				for (int i = 0; i < k; i++)
					max_size *= selected_jobs.size() - i;
				for (int i = 1; i <= k; i++)
					max_size /= i;

				std::vector<std::vector<Job_start>> result;
				result.reserve(max_size);
				std::vector<Job_start> current;
				current.reserve(k);
				combine_helper(selected_jobs, 0, k, current, result, Time_model::constants<Time>::infinity());

				assert(!result.empty());
				return result;
			}


			// assumes j is ready
			// NOTE: we don't use Interval<Time> here because the
			//       Interval c'tor sorts its arguments.
			std::pair<Time, Time> start_times(
					const State &s, const Job<Time> &j, Time t_wc, int batch_size, Time &t_high_time, Time &preempt_time) const {
				auto rt = ready_times(s, j);
				auto at = s.core_availability();
				Time est = std::max(rt.min(), at.min());

				DM("rt: " << rt << std::endl
						  << "at: " << at << std::endl);

				Time t_preempt = Time_model::constants<Time>::infinity();


				auto t_high = calculate_t_high(s, j, at.min(), t_wc + 1, batch_size); //next_higher_prio_job_ready(s, j, at.min());
				DM("t_high: " << t_high << std::endl);
				Time lst = std::min(t_wc,
									t_high - Time_model::constants<Time>::epsilon());

				lst = std::min(lst, std::max(rt.max(), s.core_availability().max()));

				// if there is a chance that the job can be dispatched,
				// we calculate the possible preemption point
				if (est <= lst)
					t_preempt = possible_preemption(est, lst + j.get_cost().max(), s, j);

//				if (t_preempt != Time_model::constants<Time>::infinity()) {
//					do {
//						auto poss_high = number_of_higher_priority(est, t_preempt, s, j);
//						if (poss_high >= num_cpus || t_preempt < s.core_availability(poss_high).max()) {
//							lst = std::min(lst, t_preempt - Time_model::constants<Time>::epsilon());
//							break;
//						} else {
//							// the job cannot be preempted.
//							// we have to check the next possible preemption
//							DM("Job cannot be preempted -> look into the next preemption point" << std::endl);
//							t_preempt = possible_preemption(t_preempt, t_preempt + j.get_cost().max(), s, j);
//							DM("Next t_preempt: " << t_preempt << std::endl);
//						}
//					} while (t_preempt < lst + j.get_cost().max());
//				}

				DM("t_preempt: " << t_preempt << std::endl);
				std::cout <<"Job" << j.get_id() << " est: " << est << " lst: " << lst << " t_preempt: " << t_preempt << std::endl;

				lst = std::min(lst, t_preempt - Time_model::constants<Time>::epsilon());
				preempt_time = t_preempt;
				t_high_time = t_high;

				DM("est: " << est << std::endl);
				DM("lst: " << lst << std::endl);

				return {est, lst};
			}

			typedef std::tuple<Job_index, Interval<Time>, Interval<Time>> Job_fin_rem; // (0) job index, (1) finish time interval, (2) remaining execution time
			void dispatch_batch(const State& s, const std::vector<Job_start>& selected_jobs) {
				// now we have to make a new state after dispatching the selected jobs
				Interval<Time> start_time = { 0, 0 };
				std::vector<Job_fin_rem> dispatched_jobs; // vector of tuple containing (0) job index, (1) finish time interval, (2) remaining execution time
				hash_value_t batch_key = 0;
				hash_value_t batch_pr_key = 0;

				for (auto it = selected_jobs.begin(); it != selected_jobs.end(); it++) {
					const Job<Time>& j = *(std::get<0>(*it));
					Interval<Time> st = std::get<1>(*it);
					const Time t_preempt = std::get<2>(*it);
					Interval<Time> cost =  s.job_preempted(index_of(j)) ? s.get_remaining_time(index_of(j)) : j.get_cost();
					start_time = Interval<Time>{ std::max(start_time.from(), st.from()), std::max(start_time.upto(), st.upto()) };

					Interval<Time> ftimes = st + cost;
					Interval<Time> remaining(0, 0);
					if (t_preempt <= ftimes.from()) {
						DM("[1] Dispatching segment: " << j << std::endl);
						remaining = { ftimes.from() - t_preempt, ftimes.upto() - t_preempt };
						ftimes = { t_preempt, t_preempt };
					}
					else if (ftimes.from() < t_preempt && t_preempt < ftimes.until()) {
						DM("[2] Dispatching segment: " << j << std::endl);
						remaining = { 0, ftimes.upto() - t_preempt };
						// since job possibly finished at time ftimes.from(),
						// if it is preempted the finish time of the previous segment is t_preempt
						ftimes = { t_preempt, t_preempt };
						// update the finish time
						update_finish_times(j, { ftimes.from(), t_preempt });
					}
					else {
						// dispatching the whole job
						DM("[3] Dispatching: " << j << std::endl);
						update_finish_times(j, ftimes);
					}


					dispatched_jobs.push_back(std::make_tuple(index_of(j), ftimes, remaining));
					if (s.job_preempted(index_of(j)) && remaining.max() == 0) {
						// the segment completed, so we have to remove it from the preempted jobs key
						batch_pr_key ^= j.get_key();
						// we have to also add it to the key of the completed jobs
						batch_key ^= j.get_key();
					}else if(!s.job_preempted(index_of(j)) && remaining.max() > 0)
						// this is a new preempted job, so we have to add it to the preempted jobs key
						batch_pr_key ^= j.get_key();
					else
						// a normal dispatched job without preemption
						batch_key ^= j.get_key();
				}

				// expand the graph, merging if possible
				const State& next = be_naive ?
					new_state(s, dispatched_jobs, start_time, batch_key, batch_pr_key) :
					new_or_merged_state(s, dispatched_jobs, start_time, batch_key, batch_pr_key);

				// make sure we didn't skip any jobs
				check_for_deadline_misses(s, next);

#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
				// create a vector of references to the selected jobs
				std::vector<Job_ref> selected_jobs_ref;
				std::vector<Interval<Time>> selected_jobs_finish;
				std::vector<Interval<Time>> selected_jobs_start;
				std::vector<Interval<Time>> selected_jobs_remaining;
				for (auto it = selected_jobs.begin(); it != selected_jobs.end(); it++) {
					const Job<Time>& j = *(std::get<0>(*it));
					selected_jobs_ref.push_back(&j);
					selected_jobs_start.push_back(std::get<1>(*it));
				}
				for(auto it = dispatched_jobs.begin(); it != dispatched_jobs.end(); it++) {
					selected_jobs_finish.push_back(std::get<1>(*it));
					selected_jobs_remaining.push_back(std::get<2>(*it));
				}
				edges.emplace_back(selected_jobs_ref, &s, &next, selected_jobs_start, selected_jobs_finish, false, selected_jobs_remaining);
#endif
				count_edge();
			}

			bool dispatch(const State &s, const Job<Time> &j, Interval<Time> st, Time t_preempt) {

				// compute range of possible finish times
				// by default we assume it is not preempted before
				Interval<Time> ftimes = st + j.get_cost();
				// check if it is preempted
				if(s.job_preempted(index_of(j))) {
					// if it is preempted, we have to use its remaining execution time
					ftimes = st + s.get_remaining_time(index_of(j));
				}

				DM("Assumed finish time: " << ftimes << std::endl);

				Interval<Time> remaining(0, 0);
				if (t_preempt <= ftimes.from()) {
					DM("[1] Dispatching segment: " << j << std::endl);
					remaining = {ftimes.from() - t_preempt, ftimes.upto() - t_preempt};
					ftimes = {t_preempt, t_preempt};
				} else if (ftimes.from() < t_preempt && t_preempt < ftimes.until()) {
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
					const State &next = be_naive ?
										new_state(s, index_of(j), predecessors_of(j),
												  st, ftimes, remaining, j.get_key()) :
										new_or_merged_state(s, index_of(j), predecessors_of(j),
															st, ftimes, remaining, j.get_key());
					// make sure we didn't skip any jobs
//                    check_for_deadline_misses(s, next);
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
					std::vector<Job_ref> selected_jobs;
					selected_jobs.push_back(&j);
					std::vector<Interval<Time>> finish_times;
					std::vector<Interval<Time>> start_times;
					std::vector<Interval<Time>> remaining_times;
					finish_times.push_back(ftimes);
					start_times.push_back(st);
					remaining_times.push_back(remaining);

					edges.emplace_back(selected_jobs, &s, &next, start_times, finish_times, true, remaining_times);

//					edges.emplace_back(&j, &s, &next,st , ftimes, true, remaining);
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
					std::vector<Job_ref> selected_jobs;
					selected_jobs.push_back(&j);
					std::vector<Interval<Time>> finish_times;
					std::vector<Interval<Time>> start_times;
					finish_times.push_back(ftimes);
					start_times.push_back(st);

					edges.emplace_back(selected_jobs, &s, &next, start_times, finish_times, false);
//					edges.emplace_back(&j, &s, &next, st, ftimes, false);
#endif
				}


				count_edge();

				return true;
			}

			void explore(const State &s) {
				bool found_one = false;

				DM("----" << std::endl);

				

				// (0) define the window of interest
				// earliest time a core is possibly available
				auto t_min = s.core_availability().min();
				// calculate how many jobs can run in parallel without influencing each other's dispatch time
				int batch_size = max_batch_size(s, t_min);
				// latest time by which a work-conserving scheduler
				// certainly schedules some job
				Time t_wc = calculate_twc(s, t_min, batch_size);

				DM("Checking state: " << s << std::endl);
				DM("t_min: " << t_min << std::endl
							 << "t_wc: " << t_wc << std::endl);

				// make a set of pointer to the eligible jobs
				// (0) job reference, (1) start time interval, (2) t_preempt
				std::vector<Job_start> eligible_jobs;
				DM("==== [1] ====" << std::endl);
				// (1) first check jobs that may be already pending
				for (const Job<Time> &j: jobs_by_win.lookup(t_min))
					if (j.earliest_arrival() <= t_min && ready(s, j) && !s.job_preempted(index_of(j))) {
						// check if this job has a feasible start-time interval
						Time t_high, t_preempt;
						auto _st = start_times(s, j, t_wc, batch_size, t_high, t_preempt);
						if (_st.first > _st.second)
							continue; // nope
						else {
							eligible_jobs.emplace_back(&j, _st, t_preempt);
						}
					}

				DM("==== [2] ====" << std::endl);
				// (2) check jobs that are released only later in the interval
				for (auto it = jobs_by_earliest_arrival.upper_bound(t_min);
					 it != jobs_by_earliest_arrival.end();
					 it++) {
					const Job<Time> &j = *it->second;
					DM(j << " (" << index_of(j) << ")" << std::endl);
					// stop looking once we've left the window of interest
					if (j.earliest_arrival() > t_wc)
						break;

					// Job could be not ready due to precedence constraints
					if (!ready(s, j))
						continue;

					// if job is previously preempted, we skip it here and we check it later
					//if (s.job_preempted(index_of(j))) //GN: aleady checked by ready(s,j)
					//	continue;

					// Since this job is released in the future, it better
					// be incomplete...
					assert(unfinished(s, j));

					// check if this job has a feasible start-time interval
					Time t_high, t_preempt;
					auto _st = start_times(s, j, t_wc, batch_size, t_high, t_preempt);
					if (_st.first > _st.second)
						continue; // nope
					else {
						eligible_jobs.emplace_back(&j, _st, t_preempt);
					}
				}

				DM("==== [3] ====" << std::endl);
				// (3) check jobs that are preempted (these jobs are assumed to be certainly released in the past)
				auto preempted_jobs = s.get_preempted_jobs();
				for (auto it = preempted_jobs.begin(); it != preempted_jobs.end(); it++) {
					const Job<Time> &j = jobs[std::get<0>(*it)];
					// check if this job has a feasible start-time interval
					Time t_high, t_preempt;
					auto _st = start_times(s, j, t_wc, batch_size, t_high, t_preempt);
					if (_st.first > _st.second)
						continue; // nope
					else {
						eligible_jobs.emplace_back(&j, _st, t_preempt);
					}
				}

				batch_size = std::min(batch_size, (int)eligible_jobs.size());
				if (batch_size > 1) {
					// now we have to find all legal combinations of the eligible jobs
					auto combinations = combine(eligible_jobs, batch_size);

					// print the combinations
					for (auto comb : combinations) {
						//						std::cout << "Combination: {";
						//						for (auto it = comb.begin(); it != comb.end(); it++) {
						//							const Job<Time> &j = **it;
						//							std::cout << j.get_id() << " ";
						//						}
						//						std::cout << "} ";
						dispatch_batch(s, comb);
						found_one = true;
					}
					//					std::cout << "--------------------" << std::endl;
					assert(found_one);
				}
				if (batch_size == 1) {
					// now we have a list of eligible jobs
					// we have to dispatch them
					for (auto it = eligible_jobs.begin(); it != eligible_jobs.end(); it++) {
						// get the job
						const Job<Time>& j = *(std::get<0>(*it));

						//						std::cout << "Regular dispatch: " << j.get_id() << std::endl;
						found_one |= dispatch(s, j, std::get<1>(*it), std::get<2>(*it));
					}
				}

				// check for a dead end
				if (!found_one && !all_jobs_scheduled(s))
					// out of options and we didn't schedule all jobs
					aborted = true;
			}

			// naive: no state merging
			void explore_naively() {
				be_naive = true;
				explore();
			}

			void explore() {
				make_initial_state();

				while (current_job_count < jobs.size()) {
					unsigned long n = 0;
#ifdef CONFIG_PARALLEL
					const auto &new_states_part = states_storage.back();
					unsigned int min_jobs = std::numeric_limits<unsigned int>::max();
					for (const States &new_states: new_states_part) {
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

					// allocate state space for next depth
					states_storage.emplace_back();


					check_depth_abort();
					check_cpu_timeout();
					if (aborted)
						break;

#ifdef CONFIG_PARALLEL
					// move states that do not have the minimum number of scheduled jobs to the next depth
					for (const States &new_states: new_states_part) {
						for (unsigned int i = 0; i < new_states.size(); i++) {
							auto &s = new_states[i];
							if (s.number_of_scheduled_jobs() != current_job_count) {
								// copy to next depth
								states().push_back(std::move(s));
//								insert_cache_state(&(*(--states().end())));
								merge_or_cache(&(*(--states().end())));
							}
						}
					}

					// then, explore the states with minimum scheduled jobs in parallel,
					// make a local copy of the number of explored states and initialize it to 0
					tbb::enumerable_thread_specific<unsigned long> partial_n;
					for (auto &ni: partial_n) {
						ni = 0;
					}
					parallel_for(new_states_part.range(),
								 [&](typename Split_states::const_range_type &r) {
									 for (auto it = r.begin(); it != r.end(); it++) {
										 const States &new_states = *it;
										 auto s = new_states.size();
										 tbb::parallel_for(tbb::blocked_range<size_t>(0, s),
														   [&](const tbb::blocked_range<size_t> &r) {
															   for (size_t i = r.begin(); i != r.end(); i++)
																   if (new_states[i].number_of_scheduled_jobs() ==
																	   current_job_count) {
																	   explore(new_states[i]);
																	   partial_n.local()++;
																   }
//											explore(new_states[i]);
														   });
									 }
								 });

					// get the number of states explored
					for (auto ni: partial_n) {
						n += ni;
					}
					// keep track of exploration front width
					width = std::max(width, n);

					num_states += n;

#else

					// Move the states with scheduled jobs more than the minimum to the next depth
					// and keep the rest state index in other_index for later exploration
					for (unsigned int i = 0; i < exploration_front.size(); i++) {
						State &s = exploration_front[i];
						if (s.number_of_scheduled_jobs() != current_job_count) {
							// copy to next depth
							states().push_back(std::move(s));
							cache_state(&(*(--states().end())));
//							merge_or_cache(&(*(--states().end())));
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
					}

					// then, explore the states with minimum scheduled jobs
					for (unsigned int i = 0; i < exploration_front.size(); i++) {
						State &s = exploration_front[i];
						if (s.number_of_scheduled_jobs() ==	current_job_count) {
							explore(s);
							n++;
							check_cpu_timeout();
							if (aborted)
								break;
						}
					}

					// keep track of exploration front width
					width = std::max(width, n);

					num_states += n;

#endif

					// clean up the state cache if necessary
					if (!be_naive)
						states_by_key.clear();

					// print number of states
					DM("Number of states: " << num_states << std::endl);
//					std::cout << "states: " << num_states << ", width: " << width << ", time: " << get_cpu_time() << std::endl;
//					std::cout << "states: " << num_states << ", width: " << width << std::endl;

#ifdef CONFIG_PARALLEL
					// propagate any updates to the response-time estimates
					for (auto& r : partial_rta)
						for (int i = 0; i < r.size(); ++i) {
							update_finish_times(rta, i, r[i]);
						}
#endif

#ifndef CONFIG_COLLECT_SCHEDULE_GRAPH
					// If we don't need to collect all states, we can remove
					// all those that we are done with, which saves a lot of
					// memory.
#ifdef CONFIG_PARALLEL
					parallel_for(states_storage.front().range(),
								 [](typename Split_states::range_type &r) {
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
								 [](typename Split_states::range_type &r) {
									 for (auto it = r.begin(); it != r.end(); it++)
										 it->clear();
								 });
#endif
					states_storage.pop_front();
				}
#endif


#ifdef CONFIG_PARALLEL
				for (auto &c: edge_counter)
					num_edges += c;
#endif
			}


#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
			friend std::ostream& operator<< (std::ostream& out,
			                                 const State_space<Time>& space)
			{
					std::map<const Schedule_state<Time>*, unsigned int> state_id;
					unsigned int i = 0;
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
						    << e.get_label()
						    << "\\nDL=" << e.deadline()
						    << "\\nES=" << e.earliest_start_time()
 						    << "\\nLS=" << e.latest_start_time()
						    << "\\nEF=" << e.earliest_finish_time()
						    << "\\nLF=" << e.latest_finish_time()
							<< "\\nRM=" << e.get_segment_remaining()
						    << "\"";
						if (e.deadline_miss_possible()) {
							out << ",color=Red,fontcolor=Red";
						}
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

namespace std {
	template<class Time>
	struct hash<Preemptive::Global::Schedule_state<Time>> {
		std::size_t operator()(Preemptive::Global::Schedule_state<Time> const &s) const {
			return s.get_key();
		}
	};
}


#endif
