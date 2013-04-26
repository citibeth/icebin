#include <glint2/MatrixMaker.hpp>
#include <glint2/IceSheet_L0.hpp>
#include <giss/IndexTranslator.hpp>
#include <glint2/HCIndex.hpp>

namespace glint2 {

// --------------------------------------------------------
/** Height Points to ice.
Uses: elev2, gcm->hpdefs, overlap */
void IceSheet_L0::compute_hp_to_ice()
{
	std::unique_ptr<giss::VectorSparseMatrix>
		overlap_m(get_overlap_m(Overlap::ICE));

	// Interpolate (for now) in height points but not X/Y
	_hp_to_ice = hp_interp(
		giss::BlitzSparseMatrix(*overlap_m), elev2, gcm->hpdefs);
}
// --------------------------------------------------------
/** Height Points to exch.
Uses: elev2, gcm->hpdefs, overlap */
void ExchSheet_L0::compute_hp_to_exch()
{
	// Create elev3, the elevation of each exchange grid cell.
	int n3 = exch.ndata();
	blitz::Array<double,1> elev3(n3);

	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell) {
		int i2 = cell->j;
		int i3 = cell->index;
		elev3[i3] = elev2[i2];
	}

	// Interpolate to exchange grid (in height point space, not X/Y).
	_hp_to_exch = hp_interp(
		giss::BlitzSparseMatrix(overlap13_m()), elev3, gcm->hpdefs);
}
// --------------------------------------------------------
/** Uses: elev2, gcm->hcmax, overlap */
std::unique_ptr<giss::VectorSparseMatrix> IceSheet_L0::ice_to_hc(
	giss::SparseAccumulator<int,double> &area1_m_hc)
{
	// Just area-weighted remap from ice back to hc grid
	auto g2_to_g1(grid2_to_grid1(
		giss::BlitzSparseMatrix(*overlap_m_hc), area1_m_hc));
// TODO: Fix up projection stuff, then un-comment line below.
//	proj_to_native(*gcm->grid1, proj, *g2_to_g1);

	return g2_to_g1;
}
// --------------------------------------------------------
/** Uses: elev2, gcm->hcmax, overlap */
std::unique_ptr<giss::VectorSparseMatrix> IceSheet_L0::exch_to_atm(
	giss::SparseAccumulator<int,double> &area1_m)
{
	// Just area-weighted remap from ice back to hc grid
	auto g2_to_g1(grid2_to_grid1(
		giss::BlitzSparseMatrix(*overlap13_m()), area1_m));
// TODO: Fix up projection stuff, then un-comment line below.
//	proj_to_native(*gcm->grid1, proj, *g2_to_g1);

	return g2_to_g1;
}
// --------------------------------------------------------

/**
@param area1_m IN/OUT: Area of each GCM cell covered by
	(non-masked-out) ice sheet.
@param area1_m_hc IN/OUT: Area of each GCM cell / height class coverd by
	(non-masked-out) ice sheet
    NOTE: Indexed in 1-D according to HCIndex convention [nhc * n1] */
void IceSheet_L0::accum_areas(
giss::SparseAccumulator<int,double> &area1_m)
{
printf("BEGIN accum_area(%s) %p\n", name.c_str(), overlap_m_hc.get());
	accum_per_row(*overlap13_m, area1_m);
printf("END accum_area(%s)\n", name.c_str());
}


// -------------------------------------------------------------


boost::function<void ()> IceSheet_L0::netcdf_define(NcFile &nc, std::string const &vname) const
{
	auto ret = IceSheet::netcdf_define(nc, vname);

	NcVar *info_var = nc.get_var((vname + ".info").c_str());
	info_var->add_att("parameterization", "L0");

	return ret;
}
// -------------------------------------------------------------


}	// namespace glint2
