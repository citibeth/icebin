#include <glint2/matrix_ops.hpp>

namespace glint2 {

// ==================================================================
// Height Classes

/** @param overlap [n1 x n2] sparse matrix */
std::unique_ptr<VectorSparseMatrix> height_classify(
VectorSparseMatrix const &overlap,
blitz::Array<double,1> const &elev2,
HeightClassifier &height_classifier)
{
	int n1 = overlap.nrow;
	int n2 = overlap.ncol;
	int nhc = height_classifier.nhc;

	// Check consistency on array sizes
	gassert(elev2.extent(0) == n2);

	std::unique_ptr<VectorSparseMatrix> ret(new VectorSparseMatrix(
		SparseDescr(nhc * n1, n2)));

	for (ii = overlap->begin(); ii != overlap->end(); ++ii) {
		int i1 = ii->row();
		int i2 = ii->col();
		int hc = height_classifier.get_hc(elev2(i2));
		ret->add(nc * n1 + i1, i2, ii->val());
	}

	return ret;
}

// ==================================================================
// Masking: Remove overlap entries for grid cells that are masked out

/** @param overlap [n1 x n2] sparse matrix
@param mask1 Row Mask: TRUE if we DON'T want to include it (same sense as numpy.ma)*/
@param mask1 Column Mask: TRUE if we DON'T want to include it (same sense as numpy.ma)*/
std::unique_ptr<VectorSparseMatrix> mask_out(
VectorSparseMatrix const &overlap,
blitz::Array<bool,1> const *mask1,
blitz::Array<bool,1> const *mask2)
{
	int n1 = overlap.nrow;
	int n2 = overlap.ncol;

	// Check consistency on array sizes
	gassert(mask1->extent(0) == n1);
	gassert(mask2->extent(0) == n2);

	std::unique_ptr<VectorSparseMatrix> ret(new VectorSparseMatrix(
		SparseDescr(n1, n2)));

	for (ii = overlap->begin(); ii != overlap->end(); ++ii) {
		int i1 = ii->row();
		if (mask1 && (*mask1)[i1]) continue;

		int i2 = ii->col();
		if (mask2 && (*mask2)[i2]) continue;

		ret->add(i1, i2, ii->val());
	}

	return ret;
}

// ==================================================================
// Area-Weighted Remapping

/** Upgrids from grid2 (ice grid) to grid1 (grid1-projected, or grid1hc-projected)
Transformation: [n2] --> [n1] */
std::unique_ptr<VectorSparseMatrix> grid2_to_grid1(
VectorSparseMatrix const &overlap)
{
	int n1 = overlap.nrow;
	int n2 = overlap.ncol;

	std::vector<double> area1(overlap.sum_per_row());
	std::vector<double> area2(overlap.sum_per_col());

	std::unique_ptr<VectorSparseMatrix> ret(new VectorSparseMatrix(
		SparseDescr(n1, n2)));

	for (ii = overlap->begin(); ii != overlap->end(); ++ii) {
		int i1 = ii->row();
		int i2 = ii->col();

		ret->add(i1, i2, ii->val() * area2[i2] / area1[i1]);
	}

	return ret;
}

std::unique_ptr<VectorSparseMatrix> grid1_to_grid2(
VectorSparseMatrix const &overlap)
{
	int n1 = overlap.nrow;
	int n2 = overlap.ncol;

	std::vector<double> area1(overlap.sum_per_row());
	std::vector<double> area2(overlap.sum_per_col());

	std::unique_ptr<VectorSparseMatrix> ret(new VectorSparseMatrix(
		SparseDescr(n2, n1)));

	for (ii = overlap->begin(); ii != overlap->end(); ++ii) {
		int i1 = ii->row();
		int i2 = ii->col();

		ret->add(i2, i1, ii->val() * area1[i1] / area2[i2]);
	}

	return ret;
}

// ==================================================================
// Geometric and Projection Error in Spherical/Cartesian Grids

// ---------------------------------------------------------------
/** Converts from values for projected grid1 to values for native grid1.
Diagonal matrix. */
void proj_to_native(Grid &grid1, VectorSparseMatrix &ret)
{
	// Set up the diagonal matrix
	std::vector<double> factors(grid1.size(), 0.0);
	for (auto cell = grid1.cells_begin(); cell != cells_end(); ++cell)
		factors.at(cell->index) = cell->proj_area / cell->native_area;

	// Multiply by it
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		ii->val() *= factors.at(ii->row());
}

void native_to_proj(Grid &grid1, VectorSparseMatrix &ret)
{
	// Set up the diagonal matrix
	std::vector<double> factors(grid1.size(), 0.0);
	for (auto cell = grid1.cells_begin(); cell != cells_end(); ++cell)
		factors.at(cell->index) = cell->native_area / cell->proj_area;

	// Multiply by it
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		ii->val() *= factors.at(ii->row());
}

// ---------------------------------------------------------------

/** Converts from values for projected grid1 to values for native grid1.
Diagonal matrix. */
std::unique_ptr<VectorSparseMatrix> proj_to_native(Grid &grid1)
{
	int n1 = grid1.ncells_full;

	std::unique_ptr<VectorSparseMatrix> ret(new VectorSparseMatrix(
		SparseDescr(n1, n1)));

	for (auto cell = grid1.cells_begin(); cell != cells_end(); ++cell) {
		ret->add(n1, n1, cell->proj_area / cell->native_area);
	}

	return ret;
}

/** Converts from values for projected grid1 to values for native grid1.
Diagonal matrix. */
std::unique_ptr<VectorSparseMatrix> native_to_proj(Grid &grid1)
{
	int n1 = grid1.ncells_full;

	std::unique_ptr<VectorSparseMatrix> ret(new VectorSparseMatrix(
		SparseDescr(n1, n1)));

	for (auto cell = grid1.cells_begin(); cell != cells_end(); ++cell) {
		ret->add(n1, n1, cell->native_area / cell->proj_area);
	}

	return ret;
}

// ==================================================================
// Interpolate Height Points in Height Space Only: (nhc, n1) --> (n2)

/** Gives weights for linear interpolation with a bunch of points.
If our point is off the end of the range, just continue the slope in extrapolation. */
static void linterp_1d(
	std::vector<double> const &xpoints,
	double xx,
	int *indices, double *weights)	// Size-2 arrays
{
	int n = xpoints.size();

	// This is the point ABOVE our value.
	// (i0 = i1 - 1, xpoints[i0] < xx <= xpoints[i1])
	// See: http://www.cplusplus.com/reference/algorithm/lower_bound/
	int i1 = lower_bound(xpoints.begin(), xpoints.end(), xx) - xpoints.begin();

	if (i1 <= 0) i1 = 1;
	if (i1 >= n) i1 = n-1;

	int i0 = i1-1;
	indices[0] = i0;
	indices[1] = i1;
	double ratio = (x - xpoints[i0]) / (xpoints[i1] - xpoints[i0]);
	weights[0] = (1.0 - ratio);
	weights[1] = ratio;
}



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
blitz::Array<double,1> elev2)
{
	int n1 = overlap.nrow;
	int n2 = overlap.ncol;
	int nhc = hpdefs.extent(0);

	// Check consistency on array sizes
	gassert(elev2.extent(0) == n2);
//	gassert(mask2.extent(0) == n2);

	//std::vector<double> area1(overlap.sum_per_row());
	std::vector<double> area2(overlap.sum_per_col());

	std::unique_ptr<VectorSparseMatrix> ret(new VectorSparseMatrix(
		SparseDescr(n2, n1 * nhc)));

	for (ii = overlap->begin(); ii != overlap->end(); ++ii) {
		int i1 = ii->row();
		int i2 = ii->col();

		double overlap_ratio = ii->val() / area2[i2];
		double elevation = std::max(elev2(i2), 0.0);

		// Interpolate in height points
		int ihps[2];
		double whps[2];	
		linterp_1d(hpdefs, elevation, ihps, whps);
		ret->add(i2, ihps[0]*n1 + i1, overlap_ratio * whps[0]);
		ret->add(i2, ihps[1]*n1 + i1, overlap_ratio * whps[1]);
	}


	return ret;
}

// ===================================================================
// CESM-Style Bi-linear Interpolation

static std::vector<double> boundaries_to_centers(
std::vector<double> const &boundaries)
{
	std::vector<double> centers;
	centers.reserve(boundaries.size()-1);
	for (int i=0; i<boundaries.size()-1; ++i)
		centers.push_back(.5 * (boundaries[i] + boundaries[i+1]));
	return centers;
}




/** Gives weights for linear interpolation with a bunch of points.
If our point is off the end of the range, just continue the slope in extrapolation. */
static void indices_1d(
	std::vector<double> const &xpoints,
	double xx,
	int *indices)	// Size-2 arrays
{
	int n = xpoints.size();

	// This is the point ABOVE our value.
	// (i0 = i1 - 1, xpoints[i0] < xx <= xpoints[i1])
	// See: http://www.cplusplus.com/reference/algorithm/lower_bound/
	int i1 = lower_bound(xpoints.begin(), xpoints.end(), xx) - xpoints.begin();

	if (i1 <= 0) i1 = 1;
	if (i1 >= n) i1 = n-1;

	int i0 = i1-1;
	indices[0] = i0;
	indices[1] = i1;
	double ratio = (x - xpoints[i0]) / (xpoints[i1] - xpoints[i0]);
	weights[0] = (1.0 - ratio);
	weights[1] = ratio;
}

// -----------------------------------------------------------------

struct {
	int i;
	int j;
	double weight;

	InterpWeight(_i, _j) : i(_i), j(_j), weight(1.0) {}
} InterpWeight;

void add_weights(VectorSparseMatrix &M, double factor, std::vector<InterpWeights> const &weight_vec, int ihp)
{
	n1 = grid1.ncells_full;
	for (auto ii = weight_vec.begin(); ii != weight_vec.end(); ++ii) {
		int i1 = grid1.ij_to_index(ii->i, ii->j);
		int i1h = ihp * n1 + i1;

		M.add(i2, i1h, factor * ii->weight);
	}
}

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
blitz::Array<bool,1> const &mask2)
{
	// ---------- Check Dimensions
	int n1 = grid1.ncells_full;
	int nhc = hpdefs.size();
	int n2 = elev2.extent(0);

	gassert(mask1.extent(1) == n1);
	gassert(mask2.extent(0) == n2);

	// --------- Compute Cell centers
	std::vector<double> lon1c(boundaries_to_centers(grid1.lonb));
	std::vector<double> lat1c(boundaries_to_centers(grid1.latb));
	std::vector<double> x2c(boundaries_to_centers(grid2.x_boundaries));
	std::vector<double> y2c(boundaries_to_centers(grid2.y_boundaries));

	std::vector<int> ilats, ilons;
	ilats.reserve(2);
	ilons.reserve(2);

	for (auto cell = grid2.cells_begin(); cell != grid2.cells_end(); ++cell) {
		int i2 = cell->index;

		// ---------- Project Cartesian cell center to sphere
		double lon2c_rad, lat2c_rad;
		transform(proj, llproj,
			x2c[cell->i], y2c[cell->j],
			lon2c_rad, lat2c_rad);
		double lon2c = lon2c_rad * R2D;
		double lat2c = lat2c_rad * R2D;

		// ---------- Find indices of nearest gridcells in lon direction
		// Note that indices may be out of range here (that's OK).
			// This is the point ABOVE our value.
			// (i0 = i1 - 1, xpoints[i0] < xx <= xpoints[i1])
			// See: http://www.cplusplus.com/reference/algorithm/lower_bound/
		int nlon = lon1c.size();
		int nearest_i[2];
		nearest_i[1] = std::lower_bound(lon1c.begin(), lon1c.end(), lon2c) - lon1c.begin();
		nearest_i[0] = nearest_i[1] - 1;

		// ----------- Find indices of nearest gridcells in lat direction
		nlat = lat1c.size();
		int nearest_j[2];
		nearest_j[1] = std::lower_bound(lat1c.begin(), lat1c.end(), lat2c) - lat1c.begin();
		nearest_j[0] = nearest_j[1] - 1;

		// ------------ Find expressions for the values at the four nearest gridcells.
		// TODO: Assumes even grid spacing for now.  Since this is extrapolation, that
		// guess is as good as any.
		// This will be either the value at that gridcell, or an average of neighbors.
		std::vector<InterpWeight> nearest_weights[2][2];
		for (int di=0; di<2; ++di) {
		for (int dj=0; dj<2; ++dj) {
			i = nearest_i[di];
			j = nearest_j[dj];
			if (i < 0 || i >= im || j < 0 || j >= jm || mask1xy(i,j)) {
				// This point is invalid.  Look for valid points among neighbors.
				nvalid = 0;
				for (ii = i-1; ii <= ii+1; ++i) {
				for (jj = j-1; jj <= jj+1; ++jj) {
					if (ii < 0) ii = im-1;
					else if (ii >= im) ii = 0;
					if (jj < 0) continue;
					else if (jj >= jm) continue;
					if (mask1xy(i,j)) continue;

					// Found a valid cell: average it.
					++nvalid;
					nearest[di][dj].push_back(InterpWeight(ii,jj));
				}}

				if (nvalid == 0) {
					fprintf("No valid neighbors for GCM grid cell (%d, %d).  "
						"Ice grid cell %d (%d,%d,%d) out of range",
						i, j, i2, cell->i, cell->j, cell->k);
					throw std::exception();
				}

				// Divide by nvalid
				double nvalid_inv = 1.0 / (double)nvalid;
				for (auto ii = nearest[di][dj].begin(); ii = nearest[di][dj].end(); ++ii)
					ii->weight *= nvalid_inv;
			} else {		// It's valid: just use the point
				nearest[di][dj].push_back(InterpWeight(i,j));
			}
		}}

		// ------------ Construct "fake" lon/lat positions for our
		// neighbor cells so bilinear interpolation will work smoothly.
		double nearest_lon[2];
		for (int k=0; k<2; ++k) {
			if (nearest_i[k] < 0) {
				nearest_lon[k] = lon1c[nlon-1] - 360.;
			} else if (nearest_i[k] >= nlon) {
				nearest_lon[k] = lon1c[0] + 360.;
			} else {
				nearest_lon[k] = lon1c[nearest_i[k]];
			}
		}
		double ratio_i = (lon2c - nearest_lon[0]) / (nearest_lon[1] - nearest_lon[0]);
		// ratio * x1 + (1-ratio) * x0

		double nearest_lat[2];
		for (int k=0; k<2; ++k) {
			if (nearest_j[k] < 0) {
				// project even spacing one gridcell beyond
				nearest_lat[k] = 2. * lat1c[0] - lat1c[1];
			} else if (nearest_j[k] >= nlat) {
				// project even spacing one gridcell beyond
				nearest_lat[k] = 2. * latc[nlat-1] - latc[nlat-2];
			} else {
				nearest_lat[k] = lat1c[nearest_j[k]];
			}
		}
		double ratio_j = (lat2c - nearest_lat[0]) / (nearest_lat[1] - nearest_lat[0]);

		// ------------ Interpolate in height classes
		int ihp[2];
		int whp[2];
		linterp_1d(hpdefs, elev2(i2), ihp, whp);

		// ------------ Assemble the interpolation
		std::vector<InterpWeight> weights;
		for (int k=0; k<2; ++k) {		// HP dimension
			add_weights(weights, whp[k] * (1.-ratio_i) * (1.-ratio_j), nearest_weights[0][0], ihp[k]);
			add_weights(weights, whp[k] * (1.-ratio_i) * (   ratio_j), nearest_weights[0][1], ihp[k]);
			add_weights(weights, whp[k] * (   ratio_i) * (1.-ratio_j), nearest_weights[1][0], ihp[k]);
			add_weights(weights, whp[k] * (   ratio_i) * (   ratio_j), nearest_weights[1][0], ihp[k]);
		}
	}
}


} // namespace glint2
