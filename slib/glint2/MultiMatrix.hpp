#include <giss/SparseMatrix.hpp>

namespace glint2 {

// --------------------------------------------------------------
/** Represents a matrix by a series of sub-matrices (almost but not quite block diagonal).
It then allows computation of y as follows:

[y]                       [x1]
[y] = diag(1/a) [M1 M2 M3][x2]
[ ]                       [x3]

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

	giss::MapSparseVector<int,double> _total_weights;

public:
	/** @param weights Scaling vector for this matrix */
	void add_matrix(
		std::unique_ptr<giss::VectorSparseMatrix> &&mat,
		giss::MapSparseVector<int,double> const &weights)
	{
		_matrices.push_back(std::move(mat));

		// --------- Compute: _total_weights += weights
		_total_weights.add(weights);
	}

	/** Computes y = diag(1/_total_weights) * M */
	void multiply(std::vector<blitz::Array<double,1>> const &xs,
		giss::VectorSparseVector<int,double> &y, bool clear_y = true, bool handle_nan = false);


	/** Computes y = diag(1/_total_weights) * M */
	void multiply(std::vector<blitz::Array<double,1>> const &xs,
		blitz::Array<double,1> &y, bool clear_y = true, bool handle_nan = false);


};



}
