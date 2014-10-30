#pragma once

#include <vector>
#include <blitz/array.h>
#include <giss/DynamicEnum.hpp>

namespace giss {

/** This is not officially part of the SparseMatrix.hpp framework.
Keeping it here is much simpler, even if not as general. */
class CSRMatrix : public std::vector<std::vector<std::pair<int, double>>> {
	typedef std::vector<std::vector<std::pair<int, double>>> super;
public:

	CSRMatrix(int nrow) : super(nrow) {}

	void add(int row, int col, double val)
		{ (*this)[row].push_back(std::make_pair(col, val)); }

	friend std::ostream &operator<<(std::ostream &out, CSRMatrix const &mat);
};

/** A relation such as y = Ax + b can be represented as:
y = A'x' where x' = [x 1].  HOWEVER... such a representation
does not allow for easy application, unless a bunch of 1's are appended to x.
Therefore, it makes sense to decompose it to remove the "unit" element.
By convention, the unit element is the last one in a dimension.

This class stores the (smaller) vector x, along with the unit vector b. */
class CSRAndUnits {
public:
	CSRMatrix mat;
	std::vector<double> units;

	CSRAndUnits(int nrow) : mat(nrow), units(nrow, 0.0) {}


	friend std::ostream &operator<<(std::ostream &out, CSRAndUnits const &matu);
};

// -------------------------------------------------
// -------------------------------------------------

/** Inputs received from the GCM may not be in the same units, etc. as
the inputs to be sent to the ice model.  This class creates a matrix
to translate between those two.

We wish to make transformations such as the following:
   y1 = x1 * by_dt               (where by_dt is 1/ the timestep)
   y2 = x2 + 273.15
   y3 = x3 + i4

These transformations can be represented via the matrix:
   y = M [x 1]
where M is:
	x1		x2	x3	x4	unit
	-----------------------------
y1 |by_dt	0	0	0	0
y2 |0		1	0	0	273.15
y3 |0		0	1	1	0

Note that the Matrix M involves variables that won't be bound until
later (eg, by_dt, where dt is the (variable) COUPLING timestep).
Therefore, the we compute M by multiplying a 3d-order tensor T with the
vector [by_dt 1]:
	M = T . [by_dt 1]

Once all variables in [by_dt 1] are bound to values (eg, on each
coupling timestep), the matrix M may be computed.  M may then be used
to compute output values (y) based on input values (x).


Note that y and x will typically be bundles of vectors (for example,
y1 is SMB, y2 is energy flux, etc).  This computation needs to be computed
for each element of the y's and x's.

The dimensions of the 3d-order tensor T are called OUTPUTS, INPUTS and
SCALARS.  This class assigns strings to each element of each
dimension, thereby allowing the user to set elements of the tensor by
string.  This is inefficient, but is only done at initialization time.
And it prevents errors. */

class VarTransformer {
public:

	enum {OUTPUTS, INPUTS, SCALARS, NDIM};	// Dimensions of our tensor

protected:
	blitz::Array<double, NDIM> _tensor;

	/** Name of each element in each dimension */
	DynamicEnum const *_ele_names[NDIM];

public:

	/** Define the name of each element in a dimension. */
	void set_names(int dim, DynamicEnum const *ele_names)
		{ _ele_names[dim] = ele_names; }

	/** Allocate the tensor, based on the names set above. */
	void allocate();

	// Accessor methods...
	DynamicEnum const &dimension(int idim) const
		{ return *_ele_names[idim]; }

	/** Set an element of the tensor, using name-based indexing.
	@return true if all OK, false on error. */
	bool set(std::string output, std::string input, std::string scalar, double val);

	/** Instantiates the scalars with specific values, and returns a 2nd-order
	matrix derived from the 3d-order tensor, in CSR format. */
	CSRAndUnits apply_scalars(
		std::vector<std::pair<std::string, double>> const &nvpairs);

	/** Apply this transformation to a vector (of things).  The vectors
	must both be pre-allocated. */
	template<type X>
	void apply(std::vector<X> const &inn, std::vector<X> &out)
	{

		giss::CSRAndUnits trans = vt.apply_scalars({
//			std::make_pair("by_dt", 1.0 / ((itime - api->itime_last) * api->dtsrc)),
			std::make_pair("unit", 1.0)});

		// Compute the variable transformation
		for (int xi=0; xi<dimension(OUTPUTS).size_nounit(); ++xi) {	// xi is index of output variable
			out[xi] = 0;

			// Consider each output variable separately...
			std::vector<std::pair<int, double>> const &row(trans.mat[xi]);
			for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
				int xj = xjj->first;		// Index of input variable
				double io_val = xjj->second;	// Amount to multiply it by
				
				out[xi] += inn[xj] * io_val;		// blitz++ vector operation
			}
		}
	}

	friend std::ostream &operator<<(std::ostream &out, VarTransformer const &vt);
};

std::ostream &operator<<(std::ostream &out, CSRMatrix const &mat);
std::ostream &operator<<(std::ostream &out, CSRAndUnits const &matu);

/** Print out the tensor as readable symbolic equations.
Used to check and debug. */
std::ostream &operator<<(std::ostream &out, VarTransformer const &vt);


}
