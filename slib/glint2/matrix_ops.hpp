#pragma once

#include <blitz/array.h>
#include <giss/SparseAccumulator.hpp>
#include <giss/SparseMatrix.hpp>
#include <giss/blitz.hpp>
#include <glint2/HeightClassifier.hpp>
#include <glint2/Grid.hpp>
#include <glint2/Grid_LonLat.hpp>
#include <glint2/Grid_XY.hpp>

namespace glint2 {

// ==================================================================
// Height Classes

/** @param overlap [n1 x n2] sparse matrix */
std::unique_ptr<giss::VectorSparseMatrix> height_classify(
giss::BlitzSparseMatrix const &overlap,
blitz::Array<double,1> const &elev2,
blitz::Array<double,1> const &hcmax);

// ==================================================================
// Masking: Remove overlap entries for grid cells that are masked out

/** @param overlap [n1 x n2] sparse matrix
@param mask1 Row Mask: TRUE if we DON'T want to include it (same sense as numpy.ma)
@param mask2 Column Mask: TRUE if we DON'T want to include it (same sense as numpy.ma)*/
std::unique_ptr<giss::VectorSparseMatrix> mask_out(
giss::BlitzSparseMatrix const &overlap,
blitz::Array<int,1> const *mask1,
blitz::Array<int,1> const *mask2);

// ==================================================================
// Area-Weighted Remapping

/** Upgrids from grid2 (ice grid) to grid1 (grid1-projected, or grid1hc-projected)
Transformation: [n2] --> [n1] */
std::unique_ptr<giss::VectorSparseMatrix> grid2_to_grid1(
	giss::BlitzSparseMatrix const &overlap,
	giss::SparseAccumulator<int,double> &area1);

void divide_by(giss::VectorSparseMatrix &mat,
	giss::SparseAccumulator<int,double> &area,
	giss::SparseAccumulator<int,double> &area_inv);

std::unique_ptr<giss::VectorSparseMatrix> grid1_to_grid2(
	giss::BlitzSparseMatrix const &overlap);

// ==================================================================
// Geometric and Projection Error in Spherical/Cartesian Grids


/** Converts from values for projected grid1 to values for native grid1.
Diagonal matrix.
@param proj Translate from native to projected coordinates.
@return Vector of the diagonals of the returned matrix. */
extern std::vector<double> proj_native_area_correct(Grid const &grid1, std::string const &sproj, std::string const &sdir);



// ==================================================================
// Interpolate Height Points in Height Space Only: (nhc, n1) --> (n2)

// ===================================================================
// CESM-Style Bi-linear Interpolation

static std::vector<double> boundaries_to_centers(
std::vector<double> const &boundaries);

// -----------------------------------------------------------------

#if 0
/** We only really expect this to work for Greenland.  Don't worry
about south pole in lon/lat coordinates and Antarctica.
[n2 x (nhc * n1)] sparse matrix

TODO: Compute center of mass explicitly so we don't have to rely on grid2 being Cartesian.  Not doing this right now because I don't know whether this interpolation method will survive.  See:
http://stackoverflow.com/questions/5271583/center-of-gravity-of-a-polygon
*/
std::unique_ptr<giss::VectorSparseMatrix> 
bilin_interp(
Grid_LonLat const &grid1,
Grid_XY const &grid2,
giss::Proj2 const &proj,
blitz::Array<double,1> const &hpdefs,
blitz::Array<double,1> const &elev2,
blitz::Array<bool,1> const &mask1,		// [n1] Shows where we will / will not expect landice
blitz::Array<bool,1> const &mask2);
#endif

} // namespace glint2

