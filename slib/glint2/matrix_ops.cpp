#include <cstdio>
#include <giss/Proj.hpp>
#include <glint2/matrix_ops.hpp>
#include <glint2/util.hpp>
#include <glint2/HCIndex.hpp>
#include <giss/constant.hpp>

namespace glint2 {

/** Divides mat /= area.
@param area (IN) Vector to divide by.  Value is moved out of this.
@param area_inv (OUT) 1/area */
void divide_by(giss::VectorSparseMatrix &mat,
	giss::SparseAccumulator<int,double> &area,
	giss::SparseAccumulator<int,double> &area_inv)
{
	// Compute 1 / area
	for (auto ii = area.begin(); ii != area.end(); ++ii)
		ii->second = 1.0d / ii->second;
	area_inv = std::move(area);

	// Divide by area.  (Now area is really equal to 1/area)
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		ii.val() *= area_inv[ii.row()];
}




} // namespace glint2
