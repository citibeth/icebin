#pragma once

#include <memory>
#include <Eigen/Sparse>
#include <giss/SparseMatrix.hpp>
#include <giss/Dict.hpp>

// ========================================================

namespace giss {

template<class SparseMatrixT>
std::unique_ptr<Eigen::SparseMatrix<double>> giss_to_Eigen(SparseMatrixT const &mat)
{
	typedef Eigen::SparseMatrix<double> EigenSM;
	std::unique_ptr<EigenSM> ret(new EigenSM(mat.nrow, mat.ncol));
	ret->setFromTriplets(
		giss::RerefIterator<decltype(mat.begin())>(mat.begin()),
		giss::RerefIterator<decltype(mat.end())>(mat.end()));
	return ret;
}	// namespace giss


inline std::unique_ptr<giss::VectorSparseMatrix> Eigen_to_giss(
	Eigen::SparseMatrix<double> const &mat)
{
	std::unique_ptr<giss::VectorSparseMatrix> ret(
		new giss::VectorSparseMatrix(giss::SparseDescr(
		mat.rows(), mat.cols())));

	int nz = mat.nonZeros();
	ret->reserve(nz);

	for (int k=0; k<mat.outerSize(); ++k) {
	for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
		ret->set(it.row(), it.col(), it.value(), SparseMatrix::DuplicatePolicy::ADD);
	}}

	return ret;
}

}
// ========================================================



// namespace Eigen {
// 
// template <class SparseMatrixT>
// class FullIterator {
// protected:
// 	SparseMatrixT *_mat;
// 	int _k;
// 	SparseMatrixT::InnerIterator _it;
// 	
// public:
// 	FullIterator(SparseMatrixT const &mat) : _mat(&mat), _k(0), _it(*_mat, k) {}
// 
// 	inline FullIterator& operator++() {
// 		++_it;
// 		if (!_it()) {
// 			++_k;
// 			if (_k < _mat->outerSize()) {
// 				_it = SparseMatrix::InnerIterator(*_mat, _k);
// 			}
// 		}
// 		return *this;
// 	}
// 
// 	inline operator bool() const { return _k < _mat->outerSize(); }
// 
// 	inline const Scalar& value() const { return it.value(); }
// 	inline Scalar& valueRef() { return it.valueRef(); }
// 
// 	inline Index index() const { return it.index(); }
// 	inline Index outer() const { return it.outer(); }
// 	inline Index row() const { return it.row(); }
// 	inline Index col() const { return it.col(); }
// };
// }	// Namespace Eign
// 
