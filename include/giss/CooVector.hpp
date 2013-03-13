#pragma once

#include <vector>
#include <algorithm>

namespace giss {

template<class IndexT, class ValT>
class CooVector : public std::vector<std::pair<IndexT, ValT>>
{
//	typedef std::vector<std::pair<IndexT, ValT>> super;

	void add(IndexT index, IndexT val)
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
