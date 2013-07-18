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



} // namespace glint2
