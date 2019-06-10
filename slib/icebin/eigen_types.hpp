#ifndef ICEBIN_EIGEN_TYPES_HPP
#define ICEBIN_EIGEN_TYPES_HPP

#include <array>
#include <ibmisc/zarray.hpp>
#include <spsparse/eigen.hpp>

namespace ibmisc {
namespace linear {
    class Weighted_Compressed;
}}

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
typedef spsparse::ConstUniverse<sparse_index_type, dense_index_type> ConstUniverseT;
typedef MakeDenseEigenT::EigenSparseMatrixT EigenSparseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic> EigenDenseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, 1> EigenColVectorT;
typedef Eigen::Matrix<val_type, 1, Eigen::Dynamic> EigenRowVectorT;
typedef Eigen::DiagonalMatrix<val_type, Eigen::Dynamic> EigenDiagonalMatrixT;
// -----------------------------------------

/** Decompresses a compressed matrix into an Eigen-type sparse matrix.
@param BvA The compressed matrix, to decompress (sparse indexing)
@param dims Dimensions to use / append when creating the decompressed matrix.
@return The decompressed matrix (dense indexing) */
EigenSparseMatrixT to_eigen_M(  // (generates in dense indexing)
ibmisc::ZArray<int,double,2> const &BvA,
std::array<SparseSetT *,2> dims);

}

#endif
