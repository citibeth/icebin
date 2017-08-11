#ifndef ICEBIN_EIGEN_TYPES_HPP
#define ICEBIN_EIGEN_TYPES_HPP

#include <spsparse/eigen.hpp>

namespace icebin {

// Types that will be used throughout as template arguments
typedef long sparse_index_type;
typedef int dense_index_type;
typedef double val_type;

// -----------------------------------------
typedef spsparse::MakeDenseEigen<sparse_index_type, val_type, 0, dense_index_type> MakeDenseEigenT;
template<int RANK>
    using TupleListT = MakeDenseEigenT::TupleListT<RANK>;
template<int RANK>
    using DenseArrayT = blitz::Array<val_type,RANK>;
typedef MakeDenseEigenT::SparseSetT SparseSetT;
typedef MakeDenseEigenT::EigenSparseMatrixT EigenSparseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic> EigenDenseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, 1> EigenColVectorT;
typedef Eigen::Matrix<val_type, 1, Eigen::Dynamic> EigenRowVectorT;
typedef Eigen::DiagonalMatrix<val_type, Eigen::Dynamic> EigenDiagonalMatrixT;
// -----------------------------------------

}

#endif
