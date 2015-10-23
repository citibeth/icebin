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

#include <vector>
#include <set>
#include <cstdio>
#include <string>
#include <giss/exit.hpp>

namespace giss {

/** Used to translated row and column matrix indices between spaces related by removal of dimensions. */
class IndexTranslator {
	std::string _name;	// For debugging
	std::vector<int> _a2b;
	std::vector<int> _b2a;
public:
	/** @param name Arbitrary name of the index to be translated, used for error messages. */
	IndexTranslator(std::string const &name) : _name(name) {}

	/** Set up the translation.
	Translation is done between indices in space A and space B.
	@param n Size of space A (indices run [0...n-1])
	@param used Set of indices that are used in space A.
	Indices in space B run [0...used.size()-1]. */
	void init(int size_a, std::set<int> const &used);

	/** Size of space A. */
	int na() const { return _a2b.size(); }

	/** Size of space B. */
	int nb() const { return _b2a.size(); }

	/** Convert an index from space A to B.
	@param a The source index, in space A.
	Input a must be in the range [0...na()-1], an exception will be thrown.
	@param check_result If true, then throw an exception on a negative return.
	@return The value in space B corresponding to input index a.  Or -1 if such a value does not exist. */
	int a2b(int a, bool check_result = true) const;


	/** Convert an index from space B to A.
	@param b The source index, in space B.
	Input b must be in the range [0...nb()-1], or an exception will be thrown.
	@param check_result If true, then throw an exception on a negative return.
	@return The value in space A corresponding to input index b.  Or -1 if such a value does not exist. */
	int b2a(int b, bool check_result = true) const;
};



inline int IndexTranslator::a2b(int a, bool check_result) const {
	if (a < 0 || a >= _a2b.size()) {
		fprintf(stderr, "a=%d is out of range (%d, %d)\n", a, 0, _a2b.size());
		giss::exit(1);
	}
	int b = _a2b[a];
	if (check_result && b < 0) {
		fprintf(stderr, "%s: a=%d produces invalid b=%d\n", _name.c_str(), a, b);
		giss::exit(1);
	}
	return b;
}

inline int IndexTranslator::b2a(int b, bool check_result) const {
	if (b < 0 || b >= _b2a.size()) {
		fprintf(stderr, "b=%d is out of range (%d, %d)\n", b, 0, _b2a.size());
		giss::exit(1);
	}
	int a = _b2a[b];
	if (check_result && a < 0) {
		fprintf(stderr, "%s: b=%d produces invalid a=%d\n", _name.c_str(), b, a);
		giss::exit(1);
	}
	return a;
}

}

#if 0
class RowColTranslator {
public :
	IndexTranslator row;
	IndexTranslator col;

	void init(
		int nrow, std::set<int> const &used_row,
		int ncol, std::set<int> const &used_col)
	{
		row(nrow, used_row);
		col(ncol, used_col) {}
}
// ------------------------------------------------------------
#endif
