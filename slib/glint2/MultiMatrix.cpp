#include <glint2/MultiMatrix.hpp>
#include <giss/sparse_vec_multiply.hpp>

namespace glint2 {

/** Computes y = diag(1/_total_area1_m) * M */
void MultiMatrix::multiply(std::vector<blitz::Array<double,1>> const &xs,
	giss::VectorSparseVector<int,double> &y, bool clear_y, bool handle_nan)
{
	if (clear_y) y.clear();

	// y = sum_i M_i x_i
	for (int i=0; i < _matrices.size(); ++i) {
		giss::multiply(*_matrices[i], xs[i], y, false);
	}

	// y /= area1_m
	giss::divide_by(y, _total_area1_m);

}


/** Computes y = diag(1/_total_area1_m) * M */
void MultiMatrix::multiply(std::vector<blitz::Array<double,1>> const &xs,
	blitz::Array<double,1> &y, bool clear_y, bool handle_nan)
{
	if (clear_y) y = 0;

	if (xs.size() != _matrices.size()) {
		fprintf(stderr, "The number of _matrices (%ld) and vectors (%ld) must match!\n", _matrices.size(), xs.size());
		throw std::exception();
	}

	// y = sum_i M_i x_i
	for (int i=0; i < _matrices.size(); ++i) {
		giss::multiply(*_matrices[i], xs[i], y, false, true);
	}

	// y /= area1_m
	for (auto ii = _total_area1_m.begin(); ii != _total_area1_m.end(); ++ii) {
		y[ii->first] /= ii->second;
	}

}

}
