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
	auto operator->() -> decltype(&*(this->IterT::operator*()))
		{ return &*(this->IterT::operator*()); }

};



/** Opposite of DerefIterator.  Requires one more
dereference than the wrapped iterator. */
template<class IterT>
class RerefIterator : public IterT {
public:
	RerefIterator(IterT const &ii) : IterT(ii) {}

	IterT &operator*() { return *this; }
	IterT *operator->() { return this; }
};



template<class superT, class KeyT, class ValT>
struct SecondIterator : public superT
{
	SecondIterator(superT const &ii) : superT(ii) {}

	ValT &operator*() const
		{ return *(superT::operator*().second); }
	ValT *operator->() const
		{ return &*(superT::operator*().second); }

	KeyT const &key() const
		{ return superT::operator*().first; }
};

template<
	template<class KeyTT, class ValTT> class BaseTpl,
	class KeyT, class ValT>
struct Dict : public BaseTpl<KeyT, std::unique_ptr<ValT>>
{
	typedef BaseTpl<KeyT, std::unique_ptr<ValT>> super;

	typedef SecondIterator<typename super::iterator, KeyT, ValT> ValIterator;
	typedef ValIterator iterator;
	ValIterator begin() { return ValIterator(super::begin()); }
	ValIterator end() { return ValIterator(super::end()); }

	typedef SecondIterator<typename super::const_iterator, KeyT, ValT> const_ValIterator;
	typedef const_ValIterator const_iterator;
	const_ValIterator begin() const { return const_ValIterator(super::begin()); }
	const_ValIterator end() const { return const_ValIterator(super::end()); }



	ValT *operator[](KeyT const &key) {
		auto ii = super::find(key);
		if (ii == super::end()) return 0;
		return &*(ii->second);
	}

	ValT const *operator[](KeyT const &key) const {
		auto ii = super::find(key);
		if (ii == super::end()) return 0;
		return &*(ii->second);
	}

	std::pair<ValT *, bool> insert(KeyT const &key, std::unique_ptr<ValT> &&valp) {
		auto ret = super::insert(std::make_pair(key, std::move(valp)));
		typename super::iterator nw_it;
			nw_it = ret.first;
		ValT *nw_valp = &*(nw_it->second);
		bool &is_inserted(ret.second);

		return std::make_pair(nw_valp, is_inserted);
	}

	std::pair<ValT *, bool> insert(KeyT const &key, ValT &&val) {
		return insert(key,
			std::unique_ptr<ValT>(new ValT(std::move(val))));
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


};	// class Dict

// -----------------------------------------
template<
	template<class KeyTT, class ValTT> class BaseTpl,
	class KeyT, class ValT>
struct Dict_Create : public Dict<BaseTpl, KeyT, ValT>
{
	typedef Dict<BaseTpl, KeyT, ValT> super;
//	typedef BaseTpl<KeyT, std::unique_ptr<ValT>> super2;

	/** Creates elements as you insert them! */
	ValT *operator[](KeyT const &key) {
		auto ii = super::find(key);
		if (ii == super::end()) {
//			std::unique_ptr<ValT> ptr(new ValT());
			auto ret = super::super::insert(std::make_pair(key, // std::move(ptr)));
				std::unique_ptr<ValT>(new ValT())));
			// typename super::super::iterator nw_it = ret.first;
			return &*(ret.first->second);
		}
		return &*(ii->second);
	}
};
// -----------------------------------------
template<class KeyT, class ValT>
class _MapDict_core : public std::map<KeyT, ValT> {};

template<class KeyT, class ValT>
class MapDict : public Dict<_MapDict_core, KeyT, ValT> {};

template<class KeyT, class ValT>
class MapDict_Create : public Dict_Create<_MapDict_core, KeyT, ValT> {};
// -----------------------------------------
template<class KeyT, class ValT>
class _HashDict_core : public std::map<KeyT, ValT> {};

template<class KeyT, class ValT>
class HashDict : public Dict<_HashDict_core, KeyT, ValT> {};

template<class KeyT, class ValT>
class HashDict_Create : public Dict_Create<_HashDict_core, KeyT, ValT> {};
// -----------------------------------------




// -----------------------------------------



}
