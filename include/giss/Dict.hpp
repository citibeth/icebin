#pragma once

#include <map>

namespace giss {


/** An iterator that derferences its result one more time.
Useful to iterate through colletions of unique_ptr as if
they're plain values. */
template<class IterT>
struct DerefIterator<IterT> : public IterT {
	DerefIterator(IterT::iterator const &ii) : IterT(ii) {}

	Vertex_L1 &operator*() { return *IterT::operator*(); }
	Vertex_L1 *operator->() { return &*IterT::operator*(); }

};


template<class KeyT, class ValT>
struct Dict<KeyT, ValT> : std::unorderd_map<KeyT, std::unique_ptr<ValT>>
{
	typedef std::unorderd_map<KeyT, std::unique_ptr<ValT>> super;

	struct ValIterator : public super::iterator
	{
		ValIterator(super::iterator const &ii) : super::iterator(ii) {}

		Vertex_L1 &operator*() { return *(super::operator*().second); }
		Vertex_L1 *operator->() { return &*(super::operator*().second); }
	};


	ValT *operator[](KeyT cont &key) {
		auto ii = find(key);
		if (ii == end()) return 0;
		return &*(ii->second);
	}

	std::pair<ValT *, bool> insert(KeyT const &key, ValT &&val) {
		std::unique_ptr<ValT> valp(new ValT(std::move(val)));
		ret = insert(std::make_pair(key, std::move(valp)));
		return std::make_pair<Vector *,bool>(&**(ret.first), ret.second);
	}

private :
	struct CmpPointers {
		bool operator()(ValT *a, ValT *b) { return *a < *b; }
	};

public :

	/** @return A sorted set of the elements. */
	std::vector<ValT *> sorted()
	{
		// Make a vector of pointers
		std::vector<ValT *> ret;
		for (auto ii = begin(); ii != end(); ++ii)
			ret.push_back(&*(ii->second));

		// Sort the vector		
		std::sort(ret.begin(), ret.end(), CmpPointers());

		return ret;
	}	


};


}
