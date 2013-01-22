#pragma once

#include <map>
#include <unordered_map>
#include <memory>
#include <string>
#include <algorithm>

namespace giss {


/** An iterator that derferences its result one more time.
Useful to iterate through colletions of unique_ptr as if
they're plain values. */
template<class IterT>
class DerefIterator : public IterT {
public:
	DerefIterator(IterT const &ii) : IterT(ii) {}

	auto operator*() -> decltype(*(this->IterT::operator*()))
		{ return *(this->IterT::operator*()); }
//	auto operator->() -> decltype(*(this->IterT::operator->()))
//		{ return &*(this->IterT::operator*()); }
	auto operator->() -> decltype(&*(this->IterT::operator*()))
		{ return &*(this->IterT::operator*()); }

};



/** An iterator that derferences its result one more time.
Useful to iterate through colletions of unique_ptr as if
they're plain values. */
template<class IterT>
class RerefIterator : public IterT {
public:
	RerefIterator(IterT const &ii) : IterT(ii) {}

	IterT &operator*() { return *this; }
	IterT *operator->() { return this; }
};



// template<class ValT>
// class DerefCmp {
// 	bool operator<(ValT const *rhs) const {
// 		return (*rhs) < (*this);
// 	}
// }

template<class KeyT, class ValT>
struct Dict : public std::unordered_map<KeyT, std::unique_ptr<ValT>>
{
	typedef std::unordered_map<KeyT, std::unique_ptr<ValT>> super;

	struct ValIterator : public super::iterator
	{
		ValIterator(typename super::iterator const &ii) : super::iterator(ii) {}

		ValT &operator*() const
			{ return *(super::iterator::operator*().second); }
		ValT *operator->() const
			{ return &*(super::iterator::operator*().second); }
#if 0
		auto operator*() -> decltype(*(super::operator*().second))
			{ return *(super::operator*().second); }
		auto operator->() -> decltype(&*(super::operator*().second))
			{ return &*(super::operator*().second); }
#endif

	};

	ValIterator begin() { return ValIterator(super::begin()); }
	ValIterator end() { return ValIterator(super::end()); }

// 	std::vector<ValT *> as_vector() {
// 		std::vector<ValT *> ret;
// 		ret.reserve(size());
// 		for (auto ii = begin(); ii != end(); ++ii) ret.push_back(&*ii);
// 		return ret;
// 	}


	struct const_ValIterator : public super::const_iterator
	{
		const_ValIterator(typename super::const_iterator const &ii) : super::const_iterator(ii) {}

		ValT const &operator*() const
			{ return *(super::const_iterator::operator*().second); }
		ValT const *operator->() const
			{ return &*(super::const_iterator::operator*().second); }
	};

	const_ValIterator begin() const { return const_ValIterator(super::begin()); }
	const_ValIterator end() const { return const_ValIterator(super::end()); }









	ValT *operator[](KeyT const &key) {
		auto ii = super::find(key);
		if (ii == super::end()) return 0;
		return &*(ii->second);
	}

	std::pair<ValT *, bool> insert(KeyT const &key, ValT &&val) {
		std::unique_ptr<ValT> valp(new ValT(std::move(val)));
		auto ret = super::insert(std::make_pair(key, std::move(valp)));
		typename super::iterator nw_it;
			nw_it = ret.first;
		ValT *nw_valp = &*(nw_it->second);
		bool &is_inserted(ret.second);

		return std::make_pair(nw_valp, is_inserted);
//		return std::make_pair<ValT *,bool>((&**(ret.first), ret.second);
	}

private :
	struct CmpPointers {
		bool operator()(ValT const *a, ValT const *b) { return *a < *b; }
	};

public :

	/** @return A sorted set of the elements. */
	std::vector<ValT *> sorted() const
	{
		// Make a vector of pointers
		std::vector<ValT *> ret;
		for (auto ii = super::begin(); ii != super::end(); ++ii)
			ret.push_back(&*(ii->second));

		// Sort the vector		
		std::sort(ret.begin(), ret.end(), CmpPointers());

		return ret;
	}	


};


}
