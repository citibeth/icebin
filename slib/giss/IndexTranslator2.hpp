#pragma once

#include <vector>
#include <set>
#include <cstdio>
#include <map>
#include <unordered_map>
#include <giss/Hash.hpp>

namespace giss {

/** Used to translated row and column matrix indices between spaces related by removal of dimensions. */
class IndexTranslator2 {

	std::map<int, size_t> _size_a;
	std::unordered_map<std::pair<int,int>, int, giss::HashPair<int,int>> _a2b;
	std::vector<std::pair<int,int>> _b2a;
public:
	/** Set up the translation.
	Translation is done between indices in space A and space B.
	@param n Size of space A (indices run [0...n-1])
	@param used Set of indices that are used in space A.
	Indices in space B run [0...used.size()-1]. */
	void init(
		std::map<int, size_t> &&size_a,
		std::set<std::pair<int,int>> const &used_a);

	size_t nindex() const { return _size_a.size(); }

	/** @return Size of space A. */
	size_t na(int index) const { return _size_a.at(index); }

	/** @return Size of space B. */
	size_t nb() const { return _b2a.size(); }

	/** Convert an index from space A to B.
	@param a The source index, in space A.
	Input a must be in the range [0...na()-1], an exception will be thrown.
	@param check_result If true, then throw an exception on a negative return.
	@return The value in space B corresponding to input index a.  Or -1 if such a value does not exist. */
	int a2b(std::pair<int,int> a, bool check_result = true) const;

	/** Convert an index from space B to A.
	@param b The source index, in space B.
	Input b must be in the range [0...nb()-1], or an exception will be thrown.
	@param check_result If true, then throw an exception on a negative return.
	@return The value in space A corresponding to input index b.  Or -1 if such a value does not exist. */
	std::pair<int,int> b2a(int b, bool check_result = true) const;
};

}
