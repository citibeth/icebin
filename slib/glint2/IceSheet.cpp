/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include <unordered_set>
#include <glint2/MatrixMaker.hpp>
#include <glint2/IceSheet.hpp>
#include <giss/memory.hpp>
#include <giss/enum.hpp>
#include <glint2/util.hpp>

namespace glint2 {

// -----------------------------------------------------
IceSheet::IceSheet() : interp_style(InterpStyle::Z_INTERP), name("icesheet") {}

IceSheet::~IceSheet() {}

/** Number of dimensions of atmosphere vector space */
size_t IceSheet::n1() const { return gcm->n1(); }


// -----------------------------------------------------
void IceSheet::clear()
{
	grid2.reset();
	exgrid.reset();
	mask2.reset();
	// elev2.clear();	// Don't know how to do this!
}
// -----------------------------------------------------
/** Call this after you've set data members, to finish construction. */
void IceSheet::realize()
{
	if (name == "") name = grid2->name;

	// Check bounds, etc.
	if (exgrid->grid1_ncells_full != gcm->grid1->ncells_full()) {
		fprintf(stderr, "Exchange Grid for %s incompatible with GCM grid: ncells_full = %d (vs %d expected)\n",
			name.c_str(), exgrid->grid1_ncells_full, gcm->grid1->ncells_full());
		throw std::exception();
	}

	if (exgrid->grid2_ncells_full != grid2->ncells_full()) {
		fprintf(stderr, "Exchange Grid for %s incompatible with Ice grid: ncells_full = %d (vs %d expected)\n",
			name.c_str(), exgrid->grid2_ncells_full, grid2->ncells_full());
		throw std::exception();
	}

	long n2 = grid2->ndata();
	if (mask2.get() && mask2->extent(0) != n2) {
		fprintf(stderr, "Mask2 for %s has wrong size: %ld (vs %ld expected)\n",
			name.c_str(), mask2->extent(0), n2);
		throw std::exception();
	}

	if (elev2.extent(0) != n2) {
		fprintf(stderr, "Elev2 for %s has wrong size: %ld (vs %ld expected)\n",
			name.c_str(), elev2.extent(0), n2);
		throw std::exception();
	}
}
// -----------------------------------------------------
std::unique_ptr<giss::VectorSparseMatrix>
IceSheet::atm_proj_correct(ProjCorrect direction)
{
	int n1 = gcm->n1();
	std::unique_ptr<giss::VectorSparseMatrix> ret(
		new giss::VectorSparseMatrix(
		giss::SparseDescr(n1, n1)));

	giss::Proj2 proj;
	gcm->grid1->get_ll_to_xy(proj, grid2->sproj);

	for (auto cell = gcm->grid1->cells_begin(); cell != gcm->grid1->cells_end(); ++cell) {
		double native_area = cell->area;
		double proj_area = area_of_proj_polygon(*cell, proj);
		ret->add(cell->index, cell->index,
			direction == ProjCorrect::NATIVE_TO_PROJ ?
				native_area / proj_area : proj_area / native_area);
	}
	ret->sum_duplicates();

	return ret;
}
// -----------------------------------------------------
/** Corrects the area1_m result */
void IceSheet::atm_proj_correct(
	giss::SparseAccumulator<int,double> &area1_m,
	ProjCorrect direction)
{
	int n1 = gcm->n1();

	giss::Proj2 proj;
	gcm->grid1->get_ll_to_xy(proj, grid2->sproj);


	for (auto ii = area1_m.begin(); ii != area1_m.end(); ++ii) {
		int i1 = ii->first;
		double area1 = ii->second;
		glint2::Cell *cell = gcm->grid1->get_cell(i1);

		double native_area = cell->area;
		double proj_area = area_of_proj_polygon(*cell, proj);
		
		// We'll be DIVIDING by area1_m, so correct the REVERSE of above.
		ii->second *= (direction == ProjCorrect::NATIVE_TO_PROJ ?
				proj_area / native_area : native_area / proj_area);
	}
}
// -----------------------------------------------------
/** Made for binding... */
static bool in_good(std::unordered_set<int> const *set, int index_c)
{
	return (set->find(index_c) != set->end());
}

void IceSheet::filter_cells1(boost::function<bool (int)> const &include_cell1)
{
	// Remove unneeded cells from exgrid
	// Figure out which cells in grid2 to keep
	std::unordered_set<int> good_index2;
	for (auto excell = exgrid->cells_begin(); excell != exgrid->cells_end(); ++excell) {
		int index1 = excell->i;
		if (include_cell1(index1)) {
			good_index2.insert(excell->j);
		} else {
			exgrid->cells_erase(excell);
		}
	}

	// Remove unneeded cells from grid2
	grid2->filter_cells(boost::bind(&in_good, &good_index2, _1));
}
// -----------------------------------------------------

// ==============================================================
// Write out the parts that this class computed --- so we can test/check them

boost::function<void ()> IceSheet::netcdf_define(NcFile &nc, std::string const &vname) const
{
	printf("Defining ice sheet %s\n", vname.c_str());

	auto one_dim = giss::get_or_add_dim(nc, "one", 1);
	NcVar *info_var = nc.add_var((vname + ".info").c_str(), ncInt, one_dim);
	info_var->add_att("name", name.c_str());

	std::vector<boost::function<void ()>> fns;

	fns.push_back(grid2->netcdf_define(nc, vname + ".grid2"));
	fns.push_back(exgrid->netcdf_define(nc, vname + ".exgrid"));

	NcDim *n2_dim = nc.add_dim((vname + ".n2").c_str(), elev2.extent(0));
	if (mask2.get()) {
		fns.push_back(giss::netcdf_define(nc, vname + ".mask2", *mask2, {n2_dim}));
	}
	fns.push_back(giss::netcdf_define(nc, vname + ".elev2", elev2, {n2_dim}));

	return boost::bind(&giss::netcdf_write_functions, fns);
}
// -------------------------------------------------------------
void IceSheet::read_from_netcdf(NcFile &nc, std::string const &vname)
{
	clear();

	printf("IceSheet::read_from_netcdf(%s) 1\n", vname.c_str());

	NcVar *info_var = nc.get_var((vname + ".info").c_str());
	name = giss::get_att(info_var, "name")->as_string(0);
	std::string sinterp_style(giss::get_att(info_var, "interp_style")->as_string(0));
	interp_style = giss::parse_enum<InterpStyle>(sinterp_style.c_str());

	grid2.reset(read_grid(nc, vname + ".grid2").release());
	exgrid = giss::dynamic_shared_cast<ExchangeGrid,Grid>(read_grid(nc, vname + ".exgrid"));
	if (giss::get_var_safe(nc, vname + ".mask2")) {
		mask2.reset(new blitz::Array<int,1>(
		giss::read_blitz<int,1>(nc, vname + ".mask2")));
	}

	elev2.reference(giss::read_blitz<double,1>(nc, vname + ".elev2"));

	printf("IceSheet::read_from_netcdf(%s) END\n", vname.c_str());
}
// ==========================================================
// ===================================================================
// CESM-Style Bi-linear Interpolation

/** Gives weights for linear interpolation with a bunch of points.
If our point is off the end of the range, just continue the slope in extrapolation.
@param xpoints This is not blitz::Array<double,1> because Blitz++ does not (yet) implement STL-compatible iterators. */
extern void linterp_1d(
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
	double ratio = (xx - xpoints[i0]) / (xpoints[i1] - xpoints[i0]);
	weights[0] = (1.0 - ratio);
	weights[1] = ratio;
}


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
	Grid const *grid1;			// Set outside the constructor
	HCIndex const *hc_index;	// Set outside the constructor
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
	for (auto ii = weight_vec.begin(); ii != weight_vec.end(); ++ii) {
		int i1 = grid1->ij_to_index(ii->i, ii->j);
		int i1h = hc_index->ik_to_index(i1, ihp);

		M->add(i2, i1h, factor * ii->weight);
	}
}
// --------------------------------------------
/** We only really expect this to work for Greenland.  Don't worry
about south pole in lon/lat coordinates and Antarctica.
@param grid2 May be ice grid or exchange grid.  But grid2, elev2 and
       mask2 all have to exist on the same grid.  It should be easy to
       convert elev2, from ice to exchange grid --- but this hasn't
       been done yet.
@return [n2 x n3] sparse matrix */
extern std::unique_ptr<giss::VectorSparseMatrix> 
bilin_interp(
MatrixMaker *gcm,
Grid const &grid1_lonlat,
Grid const &grid2,
std::vector<double> const &hpdefs,
blitz::Array<double,1> const &elev2,
blitz::Array<int,1> const *mask1,		// [n1] Shows where we will / will not expect landice
blitz::Array<int,1> const *mask2)
{
	printf("BEGIN bilin_interp(mask1=%p, mask2=%p)\n", mask1, mask2);

	// Check types
	auto grid1p = dynamic_cast<Grid_LonLat const *>(&grid1_lonlat);
	if (!grid1p) {
		fprintf(stderr, "grid1 must be of type Grid_LonLat for BILIN_INTERP\n");
		throw std::exception();
	}
	Grid_LonLat const &grid1(*grid1p);

	// Get projection
	giss::Proj2 proj;
	grid1.get_xy_to_ll(proj, grid2.sproj);

	// ---------- Check Dimensions
	long n1 = grid1.ncells_full();
//	int nhc = hpdefs.size();
	int n2 = elev2.extent(0);

	gassert(!mask1 || mask1->extent(0) == n1);
	gassert(!mask2 || mask2->extent(0) == n2);


	// --------- Compute Atmosphere Cell centers
	std::vector<double> lon1c(grid1.lonc());
	std::vector<double> lat1c(grid1.latc());

	std::vector<int> ilats, ilons;
	ilats.reserve(2);
	ilons.reserve(2);


	IJMatrixMaker mmat(giss::SparseDescr(n2, gcm->n3()));
	mmat.grid1 = &grid1;
	mmat.hc_index = &*gcm->hc_index;
	for (int i2=0; i2<grid2.ndata(); ++i2) {
//bool dolog = (i2 == 38666 || i2 == 38967);

		if (mask2 && (*mask2)(i2)) continue;	// Ignore masked-out cells

		// ---------- Get center of this cell (or point, if we're L1 grid)
		double x2ci, y2ci;
		grid2.centroid(i2, x2ci, y2ci);

		// ---------- Project Cartesian cell center to sphere
		double lon2c, lat2c;
		proj.transform(x2ci, y2ci, lon2c, lat2c);

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
			bool mask1_index1 = (mask1 && (*mask1)(index1));	// True if index1 is masked out
			if (i < 0 || i >= nlon || j < 0 || j >= nlat || mask1_index1) {
				// This point is invalid.  Look for valid points among neighbors.
				int nvalid = 0;
				for (int ii = i-1; ii <= i+1; ++ii) {
				for (int jj = j-1; jj <= j+1; ++jj) {
					if (ii < 0) ii = nlon-1;
					else if (ii >= nlon) ii = 0;
					if (jj < 0) continue;
					else if (jj >= nlat) continue;
					if (mask1_index1) continue;

					// Found a valid cell: average it.
					++nvalid;
					nearest_weights[di][dj].push_back(InterpWeight(ii,jj));
				}}

				if (nvalid == 0) {
					fprintf(stderr, "No valid neighbors for GCM grid cell (%d, %d).  "
						"Ice grid cell %d out of range",
						i, j, i2);
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
		double elev2_i2 = elev2(i2);
		linterp_1d(hpdefs, elev2_i2, ihp, whp);

		// ------------ Assemble the interpolation
		int n1 = grid1.ncells_full();
		mmat.i2 = i2;
		for (int k=0; k<2; ++k) {		// HP dimension
			mmat.add_weights(whp[k] * (1.-ratio_i) * (1.-ratio_j), nearest_weights[0][0], ihp[k]);
			mmat.add_weights(whp[k] * (1.-ratio_i) * (   ratio_j), nearest_weights[0][1], ihp[k]);
			mmat.add_weights(whp[k] * (   ratio_i) * (1.-ratio_j), nearest_weights[1][0], ihp[k]);
			mmat.add_weights(whp[k] * (   ratio_i) * (   ratio_j), nearest_weights[1][1], ihp[k]);
		}
	}

	printf("END bilin_interp()\n");

	return std::move(mmat.M);
}


// =========================================================

}	// namespace glint2
