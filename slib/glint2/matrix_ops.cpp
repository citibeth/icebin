#include <cstdio>
#include <giss/Proj.hpp>
#include <glint2/matrix_ops.hpp>
#include <glint2/util.hpp>
#include <glint2/HCIndex.hpp>
#include <giss/constant.hpp>

namespace glint2 {

// ==================================================================
// Height Classes

/** @param overlap [n1 x n2] sparse matrix */
std::unique_ptr<giss::VectorSparseMatrix> height_classify(
giss::BlitzSparseMatrix const &overlap,
blitz::Array<double,1> const &elev2,
blitz::Array<double,1> const &hcmax)
{
	HeightClassifier height_classifier(&hcmax);
	HCIndex hc_index(overlap.nrow);
	int n2 = overlap.ncol;
	int nhc = height_classifier.nhc();

	// Check consistency on array sizes
	gassert(elev2.extent(0) == n2);

	std::unique_ptr<giss::VectorSparseMatrix> ret(new giss::VectorSparseMatrix(
		giss::SparseDescr(nhc * hc_index.n1, n2)));

	for (auto ii = overlap.begin(); ii != overlap.end(); ++ii) {
		int i1 = ii.row();
		int i2 = ii.col();
		int hc = height_classifier(elev2(i2));
		int index = hc_index.ik_to_index(i1, hc);
		ret->add(index, i2, ii.val());
	}

	return ret;
}

// ==================================================================
// Masking: Remove overlap entries for grid cells that are masked out

/** @param overlap [n1 x n2] sparse matrix
@param mask1 Row Mask: TRUE if we DON'T want to include it (same sense as numpy.ma)
@param mask2 Column Mask: TRUE if we DON'T want to include it (same sense as numpy.ma)*/
std::unique_ptr<giss::VectorSparseMatrix> mask_out(
giss::BlitzSparseMatrix const &overlap,
blitz::Array<int,1> const *mask1,
blitz::Array<int,1> const *mask2)
{
	int n1 = overlap.nrow;
	int n2 = overlap.ncol;

	// Check consistency on array sizes
	if (mask1) gassert(mask1->extent(0) == n1);
	if (mask2) gassert(mask2->extent(0) == n2);

	std::unique_ptr<giss::VectorSparseMatrix> ret(new giss::VectorSparseMatrix(
		giss::SparseDescr(n1, n2)));

	for (auto ii = overlap.begin(); ii != overlap.end(); ++ii) {
		int i1 = ii.row();
		if (mask1 && (*mask1)(i1)) continue;

		int i2 = ii.col();
		if (mask2 && (*mask2)(i2)) continue;

		ret->add(i1, i2, ii.val());
	}

	return ret;
}

// ==================================================================
// Area-Weighted Remapping

/** Upgrids from grid2 (ice grid) to grid1 (grid1-projected, or grid1hc-projected)
Transformation: [n2] --> [n1]

Element (row, col) in output must be divided by area1[row] (see divide_rows())
*/
std::unique_ptr<giss::VectorSparseMatrix> grid2_to_grid1(
giss::BlitzSparseMatrix const &overlap,
giss::SparseAccumulator<int,double> &area1)
{
	int n1 = overlap.nrow;
	int n2 = overlap.ncol;

	accum_per_row(overlap, area1);

	std::unique_ptr<giss::VectorSparseMatrix> ret(new giss::VectorSparseMatrix(
		giss::SparseDescr(n1, n2)));

	for (auto ii = overlap.begin(); ii != overlap.end(); ++ii) {
		int i1 = ii.row();
		int i2 = ii.col();

		ret->add(i1, i2, ii.val());
	}

	return ret;
}

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

std::unique_ptr<giss::VectorSparseMatrix> grid1_to_grid2(
giss::BlitzSparseMatrix const &overlap)
{
	int n1 = overlap.nrow;
	int n2 = overlap.ncol;

	std::vector<double> area2(overlap.sum_per_col());

	std::unique_ptr<giss::VectorSparseMatrix> ret(new giss::VectorSparseMatrix(
		giss::SparseDescr(n2, n1)));

	for (auto ii = overlap.begin(); ii != overlap.end(); ++ii) {
		int i1 = ii.row();
		int i2 = ii.col();

		ret->add(i2, i1, ii.val() / area2[i2]);
	}

	return ret;
}

// ==================================================================
// Geometric and Projection Error in Spherical/Cartesian Grids

// ---------------------------------------------------------------
/** Converts from values for projected grid1 to values for native grid1.
Diagonal matrix.
@param proj proj[0] --> proj[1] converts from native to projected space.
@param sdir Should be "p2n" or "n2p" */
extern std::vector<double> proj_native_area_correct(Grid const &grid1, std::string const &sproj, std::string const &sdir)
{
	bool p2n = (sdir == "p2n");

	auto proj_area(grid1.get_proj_areas(sproj));
	auto native_area(grid1.get_native_areas());

	std::vector<double> num;
	std::vector<double> denom;

	if (p2n) {
		num = std::move(proj_area);
		denom = std::move(native_area);
	} else {
		denom = std::move(proj_area);
		num = std::move(native_area);
	}

	for (int i=0; i<num.size(); ++i) {
		num[i] = num[i] / denom[i];
	}

	return num;
}

// ---------------------------------------------------------------

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

// -----------------------------------------------------------------


// =====================================================================
struct InterpWeight {
	int i;
	int j;
	double weight;

	InterpWeight(int _i, int _j) : i(_i), j(_j), weight(1.0) {}
};

// A little helper class.
class IJMatrixMaker
{
public:
	Grid const *grid1;
	int i2;
	std::unique_ptr<giss::VectorSparseMatrix> M;

	IJMatrixMaker(giss::SparseDescr const &descr) :
		M(new giss::VectorSparseMatrix(descr)) {}

	void add_weights(
		double factor,
		std::vector<InterpWeight> const &weight_vec,
		int ihp);

};	// IJMatrixMaker

void IJMatrixMaker::add_weights(
	double factor,
	std::vector<InterpWeight> const &weight_vec,
	int ihp)
{
	HCIndex hc_index(grid1->ncells_full());
	for (auto ii = weight_vec.begin(); ii != weight_vec.end(); ++ii) {
		int i1 = grid1->ij_to_index(ii->i, ii->j);
		int i1h = hc_index.ik_to_index(i1, ihp);

		M->add(i2, i1h, factor * ii->weight);
	}
}
// =====================================================================
#if 0
/** We only really expect this to work for Greenland.  Don't worry
about south pole in lon/lat coordinates and Antarctica.
[n2 x (nhc * n1)] sparse matrix */
std::unique_ptr<giss::VectorSparseMatrix> 
bilin_interp(
Grid_LonLat const &grid1,
Grid_XY const &grid2,
giss::Proj2 const &proj,	// Projects XY to the Sphere
std::vector<double> const &hpdefs,
blitz::Array<double,1> const &elev2,
blitz::Array<bool,1> const &mask1,		// [n1] Shows where we will / will not expect landice
blitz::Array<bool,1> const &mask2)
{
	// ---------- Check Dimensions
	long n1 = grid1.ncells_full();
	int nhc = hpdefs.size();
	int n2 = elev2.extent(0);

	gassert(mask1.extent(1) == n1);
	gassert(mask2.extent(0) == n2);

	// --------- Compute Cell centers
	std::vector<double> lon1c(boundaries_to_centers(grid1.lonb));
	std::vector<double> lat1c(boundaries_to_centers(grid1.latb));
	std::vector<double> x2c(boundaries_to_centers(grid2.xb));
	std::vector<double> y2c(boundaries_to_centers(grid2.yb));

	std::vector<int> ilats, ilons;
	ilats.reserve(2);
	ilons.reserve(2);

	IJMatrixMaker mmat(giss::SparseDescr(n2, nhc*n1));
	mmat.grid1 = &grid1;
	for (auto cell = grid2.cells_begin(); cell != grid2.cells_end(); ++cell) {
		int i2 = cell->index;

		// ---------- Project Cartesian cell center to sphere
		double lon2c_rad, lat2c_rad;
		proj.xy2ll(
			x2c[cell->i], y2c[cell->j],
			lon2c_rad, lat2c_rad);
		double lon2c = lon2c_rad * giss::R2D;
		double lat2c = lat2c_rad * giss::R2D;

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
		int nlat = lat1c.size();
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
			int i = nearest_i[di];
			int j = nearest_j[dj];
			int index1 = grid1.ij_to_index(i, j);
			if (i < 0 || i >= nlon || j < 0 || j >= nlat || mask1(index1)) {
				// This point is invalid.  Look for valid points among neighbors.
				int nvalid = 0;
				for (int ii = i-1; ii <= i+1; ++ii) {
				for (int jj = j-1; jj <= j+1; ++jj) {
					if (ii < 0) ii = nlon-1;
					else if (ii >= nlon) ii = 0;
					if (jj < 0) continue;
					else if (jj >= nlat) continue;
					if (mask1(index1)) continue;

					// Found a valid cell: average it.
					++nvalid;
					nearest_weights[di][dj].push_back(InterpWeight(ii,jj));
				}}

				if (nvalid == 0) {
					fprintf(stderr, "No valid neighbors for GCM grid cell (%d, %d).  "
						"Ice grid cell %d (%d,%d,%d) out of range",
						i, j, i2, cell->i, cell->j, cell->k);
					throw std::exception();
				}

				// Divide by nvalid
				double nvalid_inv = 1.0 / (double)nvalid;
				for (auto ii = nearest_weights[di][dj].begin(); ii != nearest_weights[di][dj].end(); ++ii)
					ii->weight *= nvalid_inv;
			} else {		// It's valid: just use the point
				nearest_weights[di][dj].push_back(InterpWeight(i,j));
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
				nearest_lat[k] = 2. * lat1c[nlat-1] - lat1c[nlat-2];
			} else {
				nearest_lat[k] = lat1c[nearest_j[k]];
			}
		}
		double ratio_j = (lat2c - nearest_lat[0]) / (nearest_lat[1] - nearest_lat[0]);

		// ------------ Interpolate in height classes
		int ihp[2];
		double whp[2];
		linterp_1d(hpdefs, elev2(i2), ihp, whp);

		// ------------ Assemble the interpolation
		int n1 = grid1.ncells_full();
		mmat.i2 = i2;
		for (int k=0; k<2; ++k) {		// HP dimension
			mmat.add_weights(whp[k] * (1.-ratio_i) * (1.-ratio_j), nearest_weights[0][0], ihp[k]);
			mmat.add_weights(whp[k] * (1.-ratio_i) * (   ratio_j), nearest_weights[0][1], ihp[k]);
			mmat.add_weights(whp[k] * (   ratio_i) * (1.-ratio_j), nearest_weights[1][0], ihp[k]);
			mmat.add_weights(whp[k] * (   ratio_i) * (   ratio_j), nearest_weights[1][0], ihp[k]);
		}
	}

	return std::move(mmat.M);
}
#endif
// =====================================================================



} // namespace glint2
