#pragma once
#include <giss/eigen.hpp>

namespace glint2 {
	typedef Eigen::SparseMatrix<double, Eigen::RowMajor> EigenSparseMatrix;


	// Triplet(row, col, value)
	typedef Eigen::Triplet<double> Triplet;
}
