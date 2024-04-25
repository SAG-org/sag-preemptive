#ifndef INDEX_SET_H
#define INDEX_SET_H

namespace Preemptive {

	class Index_set {
	public:

		typedef std::vector<std::size_t> Set_type;

		// new empty job set
		Index_set() : the_set() {}

		// derive a new set by "cloning" an existing set and adding an index
		Index_set(const Index_set &from, std::size_t idx)
				: the_vector(std::max(from.the_vector.size(), idx + 1)) {
			std::copy(from.the_vector.begin(), from.the_vector.end(), the_vector.begin());
			the_vector[idx] = true;

			// keep the set sorted while adding the new index
			bool inserted = false;
			for (auto x: from.the_set) {
				if (x < idx) {
					the_set.push_back(x);
				} else if (x > idx && !inserted) {
					the_set.push_back(idx);
					the_set.push_back(x);
					inserted = true;
				} else {
					the_set.push_back(x);
				}
			}
			if (!inserted) {
				the_set.push_back(idx);
			}
		}

		bool operator==(const Index_set &other) const {
			return the_set == other.the_set;
		}

		bool operator!=(const Index_set &other) const {
			return the_set != other.the_set;
		}

		bool contains(std::size_t idx) const {
			return the_vector.size() > idx && the_vector[idx];
		}

		bool includes(std::vector<std::size_t> indices) const {
			for (auto i: indices)
				if (!contains(i))
					return false;
			return true;
		}

		bool is_subset_of(const Index_set &other) const {
			for (unsigned int i = 0; i < the_set.size(); i++)
				if (contains(i) && !other.contains(i))
					return false;
			return true;
		}

		std::size_t size() const {
			return the_set.size();
		}

		void add(std::size_t idx) {
			if (idx >= the_vector.size())
				the_vector.resize(idx + 1);
			the_vector[idx] = true;
			the_set.push_back(idx);
			std::sort(the_set.begin(), the_set.end());
		}

		friend std::ostream &operator<<(std::ostream &stream,
										const Index_set &s) {
			bool first = true;
			stream << "{";
			for (std::size_t i: s.the_set) {
				if (!first)
					stream << ", ";
				first = false;
				stream << i;
			}
			stream << "}";

			return stream;
		}

	private:

		Set_type the_set;
		std::vector<bool> the_vector;

		// no accidental copies
//			Index_set(const Index_set& origin) = delete;
	};
}

#endif
