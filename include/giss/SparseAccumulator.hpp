#pragma once

#include <unordered_map>

namespace giss {

template<class IndexT, class AccumT>
class SparseAccumulator : public std::unordered_map<IndexT, AccumT> {
	typedef std::unordered_map<IndexT, AccumT> super;

	void add(IndexT const &index, AccumT const &val) {
		auto ii(find(index));
		if (ii == super::end()) insert(std::make_pair(index, val));
		else (*ii) += val;
	}

	void set(IndexT const &index, AccumT const &val) {
		auto ii(find(index));
		if (ii == super::end()) insert(std::make_pair(index, val));
		else (*ii) = val;
	}


	AccumT &operator[](IndexT const &index) {
		auto ii(find(index));
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
