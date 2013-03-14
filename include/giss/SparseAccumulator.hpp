#pragma once

#include <unordered_map>

namespace giss {

template<class IndexT, class AccumT>
class SparseAccumulator : public std::unordered_map<IndexT, AccumT> {
public :
	typedef std::unordered_map<IndexT, AccumT> super;

	void add(IndexT const &index, AccumT const &val) {
		auto ii(super::find(index));
		if (ii == super::end()) super::insert(std::make_pair(index, val));
		else ii->second += val;
	}

	void set(IndexT const &index, AccumT const &val) {
		auto ii(super::find(index));
		if (ii == super::end()) super::insert(std::make_pair(index, val));
		else ii->second = val;
	}

	AccumT &operator[](IndexT const &index) {
		auto ii(super::find(index));
		if (ii == super::end()) throw std::exception();
		return ii->second;
	}

	AccumT const &operator[](IndexT const &index) const {
		auto ii(super::find(index));
		if (ii == super::end()) throw std::exception();
		return ii->second;
	}

};

template<class SparseMatrixT>
inline void accum_per_row(SparseMatrixT const &mat,
	SparseAccumulator<int,double> &accum)
{
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		accum.add(ii.row(), ii.val());
}

template<class SparseMatrixT>
inline void accum_per_col(SparseMatrixT const &mat,
	SparseAccumulator<int,double> &accum)
{
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		accum.add(ii.col(), ii.val());
}

}	// namespace giss
