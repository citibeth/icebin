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
};

/** A relation such as y = Ax + b can be represented as:
y = A'x' where x' = [x 1].  HOWEVER... such a representation
does not allow for easy application, unless a bunch of 1's are appended to x.
Therefore, it makes sense to decompose it to remove the "unit" element.
By conention, the unit element is the last one in a dimension.

This class stores the (smaller) vector x, along with the unit vector b. */
class CSRAndUnits {
public:
	CSRMatrix mat;
	std::vector<double> units;

	CSRAndUnits(int nrow) : mat(nrow), units(nrow, 0.0) {}
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

	/** Set an element of the tensor, using name-based indexing. */
	void set(std::string output, std::string input, std::string scalar, double val);

	/** Instantiates the scalars with specific values, and returns a 2nd-order
	matrix derived from the 3d-order tensor, in CSR format. */
	CSRAndUnits apply_scalars(
		std::vector<std::pair<std::string, double>> const &nvpairs);

	friend std::ostream &operator<<(std::ostream &out, VarTransformer const &vt);
};

/** Print out the tensor as readable symbolic equations.
Used to check and debug. */
std::ostream &operator<<(std::ostream &out, VarTransformer const &vt);




}
