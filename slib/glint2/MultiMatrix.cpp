#include <glint2/MultiMatrix.hpp>
#include <giss/sparse_vec_multiply.hpp>
#include <giss/exit.hpp>

namespace glint2 {

/** Computes y = diag(1/_total_weights) * Mx */
void MultiMatrix::multiply(std::vector<blitz::Array<double,1>> const &xs,
	giss::VectorSparseVector<int,double> &y, bool clear_y, bool handle_nan)
{
	if (clear_y) y.clear();

	// y = sum_i M_i x_i
	for (int i=0; i < _matrices.size(); ++i) {
		giss::multiply(*_matrices[i], xs[i], y, false);
	}

	// y /= weights
	giss::divide_by(y, _total_weights);

}


/** Computes y = diag(1/_total_weights) * Mx */
void MultiMatrix::multiply(std::vector<blitz::Array<double,1>> const &xs,
	blitz::Array<double,1> &y, bool clear_y, bool handle_nan)
{
	if (clear_y) y = 0;

	if (xs.size() != _matrices.size()) {
		fprintf(stderr, "The number of _matrices (%ld) and vectors (%ld) must match!\n", _matrices.size(), xs.size());
		giss::exit(1);
	}

	// y = sum_i M_i x_i
	for (int i=0; i < _matrices.size(); ++i) {
		giss::multiply(*_matrices[i], xs[i], y, false, true);
	}

	// y /= weights
	for (auto ii = _total_weights.begin(); ii != _total_weights.end(); ++ii) {
		y[ii->first] /= ii->second;
	}

}

}
