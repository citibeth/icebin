#pragma once

#include <giss/SparseMatrix.hpp>

namespace giss {

/// Computes y = A * x
template<class SparseMatrixT>
void multiply(SparseMatrixT const &mat,
double const * x, double *y, bool clear_y, bool handle_nan = false)
{
	int nx = mat.ncol;
	int ny = mat.nrow;
	if (clear_y) for (int iy = 0; iy < ny; ++iy) y[iy] = 0;
	for (auto ii = mat.begin(); ii != mat.end(); ++ii) {
		int ix = ii.col();
		int iy = ii.row();
		double val = ii.val() * x[ix];
		if (!handle_nan || !(std::isnan(val) || std::isinf(val)))
			y[iy] += val;
	}
}

/// Computes y = A^T * x
template<class SparseMatrixT>
void multiplyT(SparseMatrixT const &mat,
double const * x, double *y, bool clear_y, bool handle_nan = false)
{
	int nx = mat.nrow;
	int ny = mat.ncol;
	if (clear_y) for (int iy = 0; iy < ny; ++iy) y[iy] = 0;
	for (auto ii = mat.begin(); ii != mat.end(); ++ii) {
		int iy = ii.col();
		int ix = ii.row();
		double val = ii.val() * x[ix];
		if (!handle_nan || !(std::isnan(val) || std::isinf(val)))
			y[iy] += val;
	}
}

// ------------------------------------------------------------
/// Computes y = A * x
template<class SparseMatrixT>
void multiply(
	SparseMatrixT const &mat,
	blitz::Array<double,1> const &x,
	blitz::Array<double,1> &y, bool clear_y, bool handle_nan = false)
{
	int nx = mat.ncol;
	int ny = mat.nrow;
	if (clear_y) y = 0;
	for (auto ii = mat.begin(); ii != mat.end(); ++ii) {
		int ix = ii.col();
		int iy = ii.row();
		double val = ii.val() * x(ix);
		if (!handle_nan || !(std::isnan(val) || std::isinf(val)))
			y(iy) += val;
	}
}

/// Computes y = A^T * x
template<class SparseMatrixT>
void multiplyT(
	SparseMatrixT const &mat,
	blitz::Array<double,1> const &x,
	blitz::Array<double,1> &y, bool clear_y, bool handle_nan = false)
{
	int nx = mat.nrow;
	int ny = mat.ncol;
	if (clear_y) y = 0;
	for (auto ii = mat.begin(); ii != mat.end(); ++ii) {
		int iy = ii.col();
		int ix = ii.row();
		double val = ii.val() * x(ix);
		if (!handle_nan || !(std::isnan(val) || std::isinf(val)))
			y(iy) += val;
	}
}

// ------------------------------------------------------------
/// Computes y = A * x
template<class SparseMatrixT, class SparseVectorT>
void multiply(
	SparseMatrixT const &matrix,
	blitz::Array<double,1> const &x,
	SparseVectorT &y, bool clear_y, bool handle_nan = false)
{
	int nx = matrix.ncol;
	int ny = matrix.nrow;
	if (clear_y) y.clear();
	for (auto ii = matrix.begin(); ii != matrix.end(); ++ii) {
		int ix = ii.col();
		int iy = ii.row();
		double val = ii.val() * x(ix);
		if (!handle_nan || !(std::isnan(val) || std::isinf(val)))
			y.add(iy, val);
	}
}

template<class SparseMatrixT, class SparseVectorT>
void multiplyT(
	SparseMatrixT const &matrix,
	blitz::Array<double,1> const &x,
	SparseVectorT &y, bool clear_y, bool handle_nan = false)
{
	int nx = matrix.ncol;
	int ny = matrix.nrow;
	if (clear_y) y.clear();
	for (auto ii = matrix.begin(); ii != matrix.end(); ++ii) {
		int iy = ii.col();
		int ix = ii.row();
		double val = ii.val() * x(ix);
		if (!handle_nan || !(std::isnan(val) || std::isinf(val)))
			y.add(iy, val);
	}
}

}
