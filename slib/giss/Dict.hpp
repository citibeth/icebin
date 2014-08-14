/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <map>
#include <unordered_map>
#include <memory>
#include <string>
#include <algorithm>

namespace giss {


/**An interator that derferences its result one more time than the
wrapped iterator.  Implemented as a wrapper.  Useful to iterate
through colletions of unique_ptr as if they're plain values.  For
example:
@code
 MyIterator ii = ...;
 DerefIterator<MyIterator> dii(ii);
 **ii == *dii
@endcode
@param IterT The type of iterator to wrap.
@see RerefIterator
*/
template<class IterT>
class DerefIterator : public IterT {
public:
	/** Construct by wrapping an existing iterator */
	DerefIterator(IterT const &ii) : IterT(ii) {}

	auto operator*() -> decltype(*(this->IterT::operator*()))
		{ return *(this->IterT::operator*()); }
	auto operator->() -> decltype(&*(this->IterT::operator*()))
		{ return &*(this->IterT::operator*()); }

};



/** Opposite of DerefIterator.  Requires one more
dereference than the wrapped iterator.  For
example:
@code
 MyIterator ii = ...;
 RerefIterator<MyIterator> dii(ii);
 *ii == **dii
@endcode
@param IterT The type of iterator to wrap.
@see DerefIterator */
template<class IterT>
class RerefIterator : public IterT {
public:
	/** Construct by wrapping an existing iterator */
	RerefIterator(IterT const &ii) : IterT(ii) {}

	IterT &operator*() { return *this; }
	IterT *operator->() { return this; }
};


/**Wraps an iterator, accessing the <i>*->second</i> element of
that iterator.  Useful for iterating over values in a std::map or
std::unordered_set as if it were a std::vector.

@param SuperT Class that will be wrapped
@param KeyT The type of <i>SuperT->first</i>
@param ValT The type of <i>SuperT->second</i>
*/
template<class superT, class KeyT, class ValT>
struct SecondIterator : public superT
{
	/** Construct by wrapping an existing iterator */
	SecondIterator(superT const &ii) : superT(ii) {}

	/** @return the <i>->second</i> item of the underlying iterator. */
	ValT &operator*() const
		{ return *(superT::operator*().second); }
	/** @return the <i>->second</i> item of the underlying iterator. */
	ValT *operator->() const
		{ return &*(superT::operator*().second); }

	/** @return the <i>->first</i> item of the underlying iterator. */
	KeyT const &key() const
		{ return superT::operator*().first; }
};

/**A <key, value> collection of unique_ptrs, where dereferencing the unique_ptrs
is swept under the hood (as compared to standard STL collections)
@param BaseTpl The template of the underlying collection class.  Could be std::map, std::unordered_set.  Example:
@code
 typedef Dict<std::map, int, MyClass> MyDict;
 MyDict dict;
 dict.insert(17, std::unique_ptr<MyClass>(new MyClass));
 for (MyDict::iterator ii = dict.begin(); ii != dict.end(); ++ii) {
     int key = ii.key():
     MyClass &item(*ii);
 }
 MyClass *p = dict[17];
@endcode
This template is normally accessed through MapDict and HashDict, not directly.
@see MapDict, HashDict */
template<
	template<class KeyTT, class ValTT> class BaseTpl,
	class KeyT, class ValT>
class Dict
{
protected:
	typedef BaseTpl<KeyT, std::unique_ptr<ValT>> BaseClass;
	BaseClass base;

public:
	size_t size() const { return base.size();}

	/** An iterator through the values of this Dict.  Can also be used
	to get keys via ValIterator.key().
	MyDiict::ValIterator ii = myDict.begin();
	 
	@see SecondIterator */
	typedef SecondIterator<typename BaseClass::iterator, KeyT, ValT> ValIterator;

	/** Alias ValIterator so this class may be used like a standard
	STL class when iterating through, completely hiding the use of
	unique_ptr. */
	typedef ValIterator iterator;
	ValIterator begin() { return ValIterator(base.begin()); }
	ValIterator end() { return ValIterator(base.end()); }

	ValIterator erase(ValIterator &ii)
		{ return ValIterator(base.erase(ii)); }
	void clear()
		{ base.clear(); }

	/** const version of ValIterator.
	@see ValIteratorIterator */
	typedef SecondIterator<typename BaseClass::const_iterator, KeyT, ValT> const_ValIterator;

	/** Alias const_ValIterator so this class may be used like a standard
	STL class when iterating through, completely hiding the use of
	unique_ptr. */
	typedef const_ValIterator const_iterator;
	const_ValIterator begin() const { return const_ValIterator(base.begin()); }
	const_ValIterator end() const { return const_ValIterator(base.end()); }


	/** Looks up an item in the Dict by key.
	@return A pointer to the item.  Returns the null pointer if the key does not exist. */
	ValT *operator[](KeyT const &key) {
		auto ii = base.find(key);
		if (ii == base.end()) return 0;
		return &*(ii->second);
	}

	/** Looks up an item in the Dict by key.
	@return A pointer to the item.  Returns the null pointer if the key does not exist. */
	ValT const *operator[](KeyT const &key) const {
		auto ii = base.find(key);
		if (ii == base.end()) return 0;
		return &*(ii->second);
	}

	/** Inserts a new <key, value> pair.
	@param valp A unique_ptr to the value to be inserted.  The Dict
	takes ownership of the value, leaving the unique_ptr null on exit.
	@return A pair consisting of a (simple) pointer to the inserted
	value, and a boolean telling whether or not it was inserted.
	Analogous to the return values of standard STL insert()
	methods. */
	std::pair<ValT *, bool> insert(KeyT const &key, std::unique_ptr<ValT> &&valp) {
		auto ret = base.insert(std::make_pair(key, std::move(valp)));
		typename BaseClass::iterator nw_it;
			nw_it = ret.first;
		ValT *nw_valp = &*(nw_it->second);
		bool &is_inserted(ret.second);

		return std::make_pair(nw_valp, is_inserted);
	}

	/** For ValT types that support move semantics... moves the value into
	a unique_ptr, and then inserts it.
	@param key The key */
	std::pair<ValT *, bool> insert(KeyT const &key, ValT &&val) {
		return insert(key,
			std::unique_ptr<ValT>(new ValT(std::move(val))));
	}

	/** Create a new item with no-arg constructor and insert it.
	Convenience method. */
	std::pair<ValT *, bool> create(KeyT const &key) {
		return insert(key,
			std::unique_ptr<ValT>(new ValT()));
	}

private :
	struct CmpPointers {
		bool operator()(ValT const *a, ValT const *b) { return *a < *b; }
	};

public :

	/** @return A sorted vector of (simple pointers to) the values stored in the Dict. */
	std::vector<ValT *> sorted() const
	{
		// Make a vector of pointers
		std::vector<ValT *> ret;
		for (auto ii = base.begin(); ii != base.end(); ++ii)
			ret.push_back(&*(ii->second));

		// Sort the vector		
		std::sort(ret.begin(), ret.end(), CmpPointers());

		return ret;
	}	


};	// class Dict

// -----------------------------------------
/**Subclass of Dict that creates new elements via new if one attempts to
access a new key via array indexing [].
@see Dict */
template<
	template<class KeyTT, class ValTT> class BaseTpl,
	class KeyT, class ValT>
struct Dict_Create : public Dict<BaseTpl, KeyT, ValT>
{
	typedef Dict<BaseTpl, KeyT, ValT> super;

	/** Access the value of an existing key.  Creates a new item if it
	does not already exist. */
	ValT *operator[](KeyT const &key) {
		auto ii = super::base.find(key);
		if (ii == super::base.end()) {
//			std::unique_ptr<ValT> ptr(new ValT());
			auto ret = super::base.insert(std::make_pair(key, // std::move(ptr)));
				std::unique_ptr<ValT>(new ValT())));
			// typename super::base.iterator nw_it = ret.first;
			return &*(ret.first->second);
		}
		return &*(ii->second);
	}
};
// -----------------------------------------
/**(internal use only) */
template<class KeyT, class ValT>
class _MapDict_core : public std::map<KeyT, ValT> {};

/**Equal to Dict<std::map, KeyT, ValT>
@see Dict */
template<class KeyT, class ValT>
class MapDict : public Dict<_MapDict_core, KeyT, ValT> {};

/**Equal to Dict_Create<std::map, KeyT, ValT>
@see Dict_Create */
template<class KeyT, class ValT>
class MapDict_Create : public Dict_Create<_MapDict_core, KeyT, ValT> {};
// -----------------------------------------

/**(internal use only).  In spite of the name, same as _MapDict_core.
TODO: Subclass from std::unordered_map, as the nmae implies. */
template<class KeyT, class ValT>
class _HashDict_core : public std::map<KeyT, ValT> {};

/**In spite of the name, same as MapDict
TODO: Subclass from std::unordered_map, as the nmae implies. */
template<class KeyT, class ValT>
class HashDict : public Dict<_HashDict_core, KeyT, ValT> {};

/**In spite of the name, same as MapDict_Create
TODO: Subclass from std::unordered_map, as the nmae implies. */
template<class KeyT, class ValT>
class HashDict_Create : public Dict_Create<_HashDict_core, KeyT, ValT> {};
// -----------------------------------------




// -----------------------------------------



}
