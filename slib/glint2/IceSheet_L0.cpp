#include <glint2/MatrixMaker.hpp>
#include <glint2/IceSheet_L0.hpp>
#include <giss/IndexTranslator.hpp>
#include <glint2/HCIndex.hpp>

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
/** Gives weights for linear interpolation with a bunch of points.
If our point is off the end of the range, just continue the slope in extrapolation.
@param xpoints This is not blitz::Array<double,1> because Blitz++ does not (yet) implement STL-compatible iterators. */
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
	double ratio = (xx - xpoints[i0]) / (xpoints[i1] - xpoints[i0]);
	weights[0] = (1.0 - ratio);
	weights[1] = ratio;
}
// --------------------------------------------------------
// --------------------------------------------------------
/** Builds an interpolation matrix to go from height points to ice/exchange grid.
@param overlap_type Controls matrix output to ice or exchange grid. */
std::unique_ptr<giss::VectorSparseMatrix> 
IceSheet_L0::hp_interp(Overlap overlap_type)
{
	int nx = (overlap_type == Overlap::ICE ? n2() : n3());

	// Sum overlap matrix by column (ice grid cell)
	std::vector<double> area2(n2());
	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		if (masked(cell)) continue;
		area2[cell->j] += cell->area;
	}

	HCIndex hc_index(n1());
	std::unique_ptr<giss::VectorSparseMatrix> ret(new giss::VectorSparseMatrix(
		giss::SparseDescr(nx, gcm->nhp() * n1())));

	// Interpolate in the vertical
	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		if (masked(cell)) continue;

		int i1 = cell->i;
		int i2 = cell->j;
		int ix = (overlap_type == Overlap::ICE ? i2 : cell->index);

		double overlap_ratio = cell->area / area2[i2];
		double elevation = std::max(elev2(i2), 0.0);

		// Interpolate in height points
		int ihps[2];
		double whps[2];	
		linterp_1d(gcm->hpdefs, elevation, ihps, whps);
//		printf("%d %f %f (%d : %f), (%d : %f)\n", i2, elev2(i2), overlap_ratio, ihps[0], whps[0], ihps[1], whps[1]);
		ret->add(ix, hc_index.ik_to_index(i1, ihps[0]),
			overlap_ratio * whps[0]);
		ret->add(ix, hc_index.ik_to_index(i1, ihps[1]),
			overlap_ratio * whps[1]);
	}


	return ret;
}
// --------------------------------------------------------
// --------------------------------------------------------
std::unique_ptr<giss::VectorSparseMatrix> IceSheet_L0::hp_to_ice()
{
	// Interpolate (for now) in height points but not X/Y
	return hp_interp(Overlap::ICE);
}
// --------------------------------------------------------
std::unique_ptr<giss::VectorSparseMatrix> IceSheet_L0::hp_to_atm(
	giss::SparseAccumulator<int,double> &area1_m)
{
	// ============= hp_to_exch
	// Interpolate (for now) in height points but not X/Y
	auto hp_to_exch(hp_interp(Overlap::EXCH));

	// ============= exch_to_atm (with area1 scaling factor)
	// Area-weighted remapping from exchange to GCM grid is equal
	// to scaled version of overlap matrix.
	std::unique_ptr<giss::VectorSparseMatrix> exch_to_atm(
		new giss::VectorSparseMatrix(giss::SparseDescr(n1(), n3())));
	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		if (masked(cell)) continue;

		// Exchange Grid is in Cartesian coordinates
		double area = cell->area;	// Computed in ExchangeGrid::overlap_callback()
		exch_to_atm->add(cell->i, cell->index, area);
		area1_m.add(cell->i, area);
	}

// TODO: Fix up projection stuff, then un-comment line below.
//	proj_to_native(*gcm->grid1, proj, *exch_to_atm);

	// ============= hp_to_atm
	return multiply(*exch_to_atm, *hp_to_exch);
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
// ==============================================================
boost::function<void ()> IceSheet_L0::netcdf_define(NcFile &nc, std::string const &vname) const
{
	auto ret = IceSheet::netcdf_define(nc, vname);

	NcVar *info_var = nc.get_var((vname + ".info").c_str());
	info_var->add_att("parameterization", "L0");

	return ret;
}
// -------------------------------------------------------------


}	// namespace glint2
