#pragma once

#include <unordered_map>

namespace giss {

template<class IndexT, class AccumT>
class SparseAccumulator : public std::unordered_map<IndexT, AccumT> {
	void add(IndexT const &index, AccumT const &val) {
		auto ii(find(index));
		if (ii == end()) insert(std::make_pair(index, val));
		else (*ii) += val;
	}

	void set(IndexT const &index, AccumT const &val) {
		auto ii(find(index));
		if (ii == end()) insert(std::make_pair(index, val));
		else (*ii) = val;
	}


	AccumT &operator[](IndexT const &index) {
		auto ii(find(index));
		if (ii == end()) throw std::exception();
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
inline void accum_per_row(SparseMatrixT const &mat,
	SparseAccumulator<int,double> &accum)
{
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		accum.add(ii.col(), ii.val());
}

#if 0
template<class SparseMatrixT>
inline void divide_rows(SparseMatrixT &mat,
	SparseAccumulator<int,double> const &div)
{
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		ii.val() /= div[ii.row()];
}

template<class SparseMatrixT>
inline void divide_cols(SparseMatrixT &mat,
	SparseAccumulator<int,double> const &div)
{
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		ii.val() /= div[ii.col()];
}
#endif


}	// namespace giss
