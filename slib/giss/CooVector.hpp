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
#include <algorithm>

namespace giss {


/** Implements a "sparse vector", stored as an array of <index, value>
pairs.  This is analogous to a coordinate-style sparse matrix
representation (VectorSparseMatrix).  The SparseAccumulator class does
a similar job, but it is analogous to a lookup-based sparse matrix
representation (MapSparseMatrix).
@see SparseAccumulator, VectorSparseMatrix, MapSparseMatrix */
template<class IndexT, class ValT>
class CooVector : public std::vector<std::pair<IndexT, ValT>>
{
public:
//	typedef std::vector<std::pair<IndexT, ValT>> super;

	void add(IndexT const &index, ValT const &val)
		{ this->push_back(std::make_pair(index, val)); }

	void sort()
		{ std::sort(this->begin(), this->end()); }

	void sum_duplicates();
};
// ---------------------------------------------------
template<class IndexT, class ValT>
void CooVector<IndexT,ValT>::sum_duplicates()
{
	sort();

	// Scan through, overwriting our array
	// New array will never be bigger than original
	int j=-1;		// Last-written item in output array
	int last_index = -1;	// Last index we saw in input array
	for (int i=0; i<this->size(); ++i) {
		int index = (*this)[i].first;
		int val = (*this)[i].second;

		if (index == last_index) {
			(*this)[j].second += (*this)[i].second;
		} else {
			++j;
			(*this)[j].second = (*this)[i].second;
			last_index = index;
		}
	}
	int n = j+1;	// Size of output array

	this->resize(n);
}




}
