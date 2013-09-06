#include <cstdio>
#include <glint2/MatrixMaker.hpp>
#include <glint2/IceSheet_L0.hpp>
#include <giss/IndexTranslator.hpp>
#include <glint2/HCIndex.hpp>
#include <giss/enum.hpp>

namespace glint2 {

// --------------------------------------------------------
/** Tells whether a cell in the exchange grid is masked out or not */
bool IceSheet_L0::masked(giss::HashDict<int, Cell>::iterator const &it)
{
	if (gcm->mask1.get() && (*gcm->mask1)(it->i)) return true;
	if (mask2.get() && (*mask2)(it->j)) return true;
	return false;
}
/** Tells whether a cell in the exchange grid is masked out or not */
bool IceSheet_L0::masked(giss::HashDict<int, Cell>::const_iterator const &it)
{
	if (gcm->mask1.get() && (*gcm->mask1)(it->i)) return true;
	if (mask2.get() && (*mask2)(it->j)) return true;
	return false;
}
// --------------------------------------------------------

/** Does elevation-class-style interpolation on height points.  Assumes
elevation class boundaries midway between height points.
@return Index of point in xpoints[] array that is closes to xx. */
static int nearest_1d(
	std::vector<double> const &xpoints,
	double xx)
{
	int n = xpoints.size();

	// This is the point ABOVE our value.
	// (i0 = i1 - 1, xpoints[i0] < xx <= xpoints[i1])
	// See: http://www.cplusplus.com/reference/algorithm/lower_bound/
	int i1 = lower_bound(xpoints.begin(), xpoints.end(), xx) - xpoints.begin();

	// Convert to point NEAREST ot ours
	if (i1 <= 0) return 0;
	else if (i1 >= n) return n-1;
	else {
		int i0 = i1-1;

		// Distance to i0 vs i1
		double d0 = std::abs(xx - xpoints[i0]);
		double d1 = std::abs(xpoints[i1] - xx);

		if (d0 <= d1) return i0;
		return i1;
	}
}


// --------------------------------------------------------
// =========================================================
// Stuff for bilin_interp()


// --------------------------------------------------------
/** Builds an interpolation matrix to go from height points to ice/exchange grid.
@param dest Controls matrix output to ice or exchange grid. */
std::unique_ptr<giss::VectorSparseMatrix> 
IceSheet_L0::hp_to_iceexch(IceExch dest)
{
printf("BEGIN hp_interp(%)\n", dest.str());
	if (interp_style == InterpStyle::BILIN_INTERP) {

		if (dest == IceExch::EXCH) {
			fprintf(stderr, "bilin_interp() does not (yet) interpolate to exchange grid, only ice grid.  This functionality would not be hard to implement.  But it's not needed because the exchange grid is only used in cases where the RM matrix can be kept local.  bilin_interp() introduces inherent non-locality into the matrix M.\n");
			throw std::exception();
		}

		return bilin_interp(gcm, *gcm->grid1, *grid2,
//			dest == IceExch::ICE ? *grid2 : *exgrid,
			gcm->hpdefs, elev2, &*gcm->mask1, &*mask2);
	}

	// Sum overlap matrix by column (ice grid cell)
	std::vector<double> area2(n2());
	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		if (masked(cell)) continue;
		area2[cell->j] += cell->area;
	}


printf("MID hp_interp(%s)\n", dest.str());

	int nx = niceexch(dest);

	std::unique_ptr<giss::VectorSparseMatrix> ret(new giss::VectorSparseMatrix(
		giss::SparseDescr(nx, gcm->n3())));

	// Interpolate in the vertical
	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		if (masked(cell)) continue;

		int const i1 = cell->i;
		int const i2 = cell->j;
		int const i4 = cell->index;
		int ix = (dest == IceExch::ICE ? i2 : i4);

		double overlap_ratio =
			(dest == IceExch::ICE ? cell->area / area2[i2] : 1.0);
		double elevation = std::max(elev2(i2), 0.0);
//elevation = 800;

		// Interpolate in height points
		switch(interp_style.index()) {
			case InterpStyle::Z_INTERP : {
				int ihps[2];
				double whps[2];
				linterp_1d(gcm->hpdefs, elevation, ihps, whps);
//printf("LINTERP %d %f %f (%d : %f), (%d : %f)\n", i2, elev2(i2), overlap_ratio, ihps[0], whps[0], ihps[1], whps[1]);
				ret->add(ix, gcm->hc_index->ik_to_index(i1, ihps[0]),
					overlap_ratio * whps[0]);
				ret->add(ix, gcm->hc_index->ik_to_index(i1, ihps[1]),
					overlap_ratio * whps[1]);
			} break;
			case InterpStyle::ELEV_CLASS_INTERP : {
				int ihps0 = nearest_1d(gcm->hpdefs, elevation);
				ret->add(ix, gcm->hc_index->ik_to_index(i1, ihps0),
					overlap_ratio);
			} break;
		}
	}

printf("END hp_interp(%s)\n", dest.str());

	return ret;
}
// --------------------------------------------------------
// --------------------------------------------------------
blitz::Array<double,1> const IceSheet_L0::ice_to_interp(
	blitz::Array<double,1> const &f2)
{
	if (interp_grid == IceExch::ICE) return f2;

	blitz::Array<double,1> f4(n4());
	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		if (masked(cell)) continue;
		// cell->i = index in atmosphere grid
		int i2 = cell->j;		// index in ice grid
		int i4 = cell->index; 	// index in exchange grid
		f4(i4) = f2(i2);
	}
	return f4;
}
// --------------------------------------------------------
std::unique_ptr<giss::VectorSparseMatrix> IceSheet_L0::hp_to_projatm(
	giss::SparseAccumulator<int,double> &area1_m)
{
printf("BEGIN IceSheet_L0::hp_to_projatm %ld %ld\n", n1(), n4());

	// ============= hp_to_interp
	// Interpolate (for now) in height points but not X/Y
	auto hp_to_interp(hp_to_iceexch(interp_grid));

	// ============= interp_to_projatm (with area1 scaling factor)
	auto interp_to_projatm(iceexch_to_projatm(area1_m, interp_grid));

	// ============= hp_to_projatm
	auto ret(multiply(*interp_to_projatm, *hp_to_interp));
printf("END IceSheet_L0::hp_to_projatm()\n");
	return ret;
}
// --------------------------------------------------------
std::unique_ptr<giss::VectorSparseMatrix> IceSheet_L0::iceexch_to_projatm(
	giss::SparseAccumulator<int,double> &area1_m,
	IceExch src)
{
printf("BEGIN IceSheet_L0::ice_to_projatm %ld %ld\n", n1(), n4());

	// ============= ice_to_atm (with area1 scaling factor)
	// Area-weighted remapping from exchange to atmosphere grid is equal
	// to scaled version of overlap matrix.
	int nx = niceexch(src);
	std::unique_ptr<giss::VectorSparseMatrix> ice_to_projatm(
		new giss::VectorSparseMatrix(giss::SparseDescr(n1(), nx)));
	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		if (masked(cell)) continue;

		// Exchange Grid is in Cartesian coordinates
		// cell->i = index in atmosphere grid
		// cell->j = index in ice grid
		// cell->index = index in exchange grid
		double area = cell->area;	// Computed in ExchangeGrid::overlap_callback()
		ice_to_projatm->add(cell->i,
			src == IceExch::ICE ? cell->j : cell->index,
			area);
		area1_m.add(cell->i, area);
	}

	//ice_to_projatm->sum_duplicates();

	return ice_to_projatm;
}
// --------------------------------------------------------
/**
@param area1_m IN/OUT: Area of each GCM cell covered by
	(non-masked-out) ice sheet. */
void IceSheet_L0::accum_areas(
giss::SparseAccumulator<int,double> &area1_m)
{
printf("BEGIN accum_area(%s)\n", name.c_str());

	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		if (masked(cell)) continue;

		area1_m.add(cell->i, cell->area);
	}
printf("END accum_area(%s)\n", name.c_str());
}
// -------------------------------------------------------------

boost::function<void ()> IceSheet_L0::netcdf_define(NcFile &nc, std::string const &vname) const
{
	auto ret = IceSheet::netcdf_define(nc, vname);

	NcVar *info_var = nc.get_var((vname + ".info").c_str());
	info_var->add_att("parameterization", "L0");
	info_var->add_att("interp_grid", interp_grid.str());

	return ret;
}
// -------------------------------------------------------------
void IceSheet_L0::read_from_netcdf(NcFile &nc, std::string const &vname)
{
	IceSheet::read_from_netcdf(nc, vname);

	NcVar *info_var = nc.get_var((vname + ".info").c_str());
	std::string sinterp_grid(giss::get_att(info_var, "interp_grid")->as_string(0));

	interp_grid = giss::parse_enum<IceExch>(sinterp_grid.c_str());
}


}	// namespace glint2
