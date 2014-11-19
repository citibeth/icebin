#include <giss/SparseMatrix.hpp>

namespace glint2 {

// --------------------------------------------------------------
/** Represents a matrix by a series of sub-matrices (almost but not quite block diagonal).
It then allows computation of y as follows:

[  ]            [M1][x1]
[y] = diag(1/a) [M2][x2]
[ ]             [M3][x3]

where:
   M1..Mn = matrices (1 per ice sheet)
   x1..xn = vectors (1 field on each ice sheet)
       y  = result (on eg atmosphere grid)
       a  = scaling vector (on same grid as y) = sum(a1 ... an)

   a1..an = scalaing vector for each matrix (dimensions same as y)

*/
class MultiMatrix
{
	std::vector<std::unique_ptr<giss::VectorSparseMatrix>> _matrices;

	giss::MapSparseVector<int,double> _total_area1_m;

public:
	/** @param area1_m Scaling vector for this matrix */
	void add_matrix(
		std::unique_ptr<giss::VectorSparseMatrix> &&mat,
		giss::MapSparseVector<int,double> const &area1_m)
	{
		_matrices.push_back(std::move(mat));

		// --------- Compute: _total_area1_m += area1_m
		_total_area1_m.add(area1_m);
	}

#if 0
	/** Computes y = diag(1/_total_area1_m) * M */
	void multiply(std::vector<blitz::Array<double,1>> const &xs,
		blitz::Array<double,1> &y, bool clear_y = true)
	{
		if (clear_y) y = 0;

		if (xs.size() != _matrices.size()) {
			fprintf(stderr, "The number of _matrices (%ld) and vectors (%ld) must match!\n", _matrices.size(), xs.size());
			throw std::exception();
		}

		// y = sum_i M_i x_i
		for (int i=0; i < _matrices.size(); ++i) {
			giss::multiply(*_matrices[i], xs[i], y, false);
		}

		// y /= area1_m
		for (auto ii = _total_area1_m.begin(); ii != _total_area1_m.end(); ++ii) {
			y[ii->first] /= ii->second;
		}

	}
#endif

	/** Computes y = diag(1/_total_area1_m) * M */
	void multiply(std::vector<blitz::Array<double,1>> const &xs,
		giss::VectorSparseVector<int,double> &y, bool clear_y = true)
	{
		if (clear_y) y.clear();

		// y = sum_i M_i x_i
		for (int i=0; i < _matrices.size(); ++i) {
			giss::multiply(*_matrices[i], xs[i], y, false);
		}

		// y /= area1_m
		giss::divide_by(y, _total_area1_m);

	}


};



}
