#ifndef INDEX_SET_H
#define INDEX_SET_H

namespace NP {

		class Index_set
		{
			public:

			typedef std::vector<std::size_t> Set_type;

			// new empty job set
			Index_set() : the_set() {}

			// derive a new set by "cloning" an existing set and adding an index
			Index_set(const Index_set& from, std::size_t idx)
			: the_set(from.the_set)
			{
				the_set.push_back(idx);
				std::sort(the_set.begin(), the_set.end());
			}

			bool operator==(const Index_set &other) const
			{
				return the_set == other.the_set;
			}

			bool operator!=(const Index_set &other) const
			{
				return the_set != other.the_set;
			}

			bool contains(std::size_t idx) const
			{
				return std::find(the_set.begin(), the_set.end(), idx) != the_set.end();
			}

			bool includes(std::vector<std::size_t> indices) const
			{
				for (auto i : indices)
					if (!contains(i))
						return false;
				return true;
			}

			bool is_subset_of(const Index_set& other) const
			{
				for (unsigned int i = 0; i < the_set.size(); i++)
					if (contains(i) && !other.contains(i))
						return false;
				return true;
			}

			std::size_t size() const
			{
				return the_set.size();
			}

			void add(std::size_t idx)
			{
				the_set.push_back(idx);
				std::sort(the_set.begin(), the_set.end());
			}

			friend std::ostream& operator<< (std::ostream& stream,
			                                 const Index_set& s)
			{
				bool first = true;
				stream << "{";
				for (auto i = 0; i < s.the_set.size(); i++) {
					if (!first)
						stream << ", ";
					first = false;
					stream << s.the_set[i];
				}
				stream << "}";

				return stream;
			}

			private:

			Set_type the_set;

			// no accidental copies
//			Index_set(const Index_set& origin) = delete;
		};
}

#endif
