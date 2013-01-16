#pragma once

#include <glint2/eigen.hpp>

namespace glint2 {

// ==================================================================
// Height Classes

/** @param overlap [n1 x n2] sparse matrix */
std::unique_ptr<VectorSparseMatrix> height_classify(
VectorSparseMatrix const &overlap,
blitz::Array<double,1> const &elev2,
HeightClassifier &height_classifier);

// ==================================================================
// Masking: Remove overlap entries for grid cells that are masked out

/** @param overlap [n1 x n2] sparse matrix
@param mask1 Row Mask: TRUE if we DON'T want to include it (same sense as numpy.ma)*/
@param mask1 Column Mask: TRUE if we DON'T want to include it (same sense as numpy.ma)*/
std::unique_ptr<VectorSparseMatrix> mask_out(
VectorSparseMatrix const &overlap,
blitz::Array<bool,1> const *mask1,
blitz::Array<bool,1> const *mask2);

// ==================================================================
// Area-Weighted Remapping

/** Upgrids from grid2 (ice grid) to grid1 (grid1-projected, or grid1hc-projected)
Transformation: [n2] --> [n1] */
std::unique_ptr<VectorSparseMatrix> grid2_to_grid1(
VectorSparseMatrix const &overlap);

std::unique_ptr<VectorSparseMatrix> grid1_to_grid2(
VectorSparseMatrix const &overlap);

// ==================================================================
// Geometric and Projection Error in Spherical/Cartesian Grids

void proj_to_native(Grid &grid1, VectorSparseMatrix &ret);
void native_to_proj(Grid &grid1, VectorSparseMatrix &ret);



/** Converts from values for projected grid1 to values for native grid1.
Diagonal matrix. */
std::unique_ptr<VectorSparseMatrix> proj_to_native(Grid &grid1);

/** Converts from values for projected grid1 to values for native grid1.
Diagonal matrix. */
std::unique_ptr<VectorSparseMatrix> native_to_proj(Grid &grid1);

// ==================================================================
// Interpolate Height Points in Height Space Only: (nhc, n1) --> (n2)



/**
		overlap[n1, n2] (scipy.sparse.coo_matrix):
			The overlap matrix between grid1 (GCM) and grid2 (ice).
			Unused cells in grid2 (and maybe also grid1h) should already be masked out.
		_elev1h[nhc, n1] (np.array):
			Set of elevation points in each grid cell we're computing on.
			Frequently, elevation points are the same for all grid cells.
			(may be any shape, as long as shape[0]=nhc and shape[1:] = n1)
		_elev2[n2] (np.array):
			Elevation of each ice grid cell (or grid point)
			(may be any shape, as long as it has n2 elements)
		_mask2[n2] (np.array, dtype=bool):
			True for gridcells in ice grid that have ice.
			(may be any shape, as long as it has n2 elements)
			NOTE: The sense of this mask is SAME that used in numpy.ma (true = masked out)

*/
std::unique_ptr<VectorSparseMatrix> 
hc_interp(VectorSparseMatrix const &overlap,
std::vector<double> const &hpdefs,
blitz::Array<double,1> elev2);

// ===================================================================
// CESM-Style Bi-linear Interpolation

static std::vector<double> boundaries_to_centers(
std::vector<double> const &boundaries);



/** Gives weights for linear interpolation with a bunch of points.
If our point is off the end of the range, just continue the slope in extrapolation. */
static void indices_1d(
	std::vector<double> const &xpoints,
	double xx,
	int *indices);	// Size-2 arrays

// -----------------------------------------------------------------

/** We only really expect this to work for Greenland.  Don't worry
about south pole in lon/lat coordinates and Antarctica.
[n2 x (nhc * n1)] sparse matrix */
bilin_interp(
Grid_LonLat const &grid1,
Grid_XY const &grid2,
Proj const &proj,
std::vector<double> const &hpdefs,
blitz::Array<double,1> const &elev2,
blitz::Array<bool,1> const &mask1,		// [n1] Shows where we will / will not expect landice
blitz::Array<bool,1> const &mask2);

} // namespace glint2

