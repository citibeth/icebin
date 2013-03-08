#include <glint2/MatrixMaker.hpp>
#include <glint2/IceSheet_L0.hpp>
#include <giss/IndexTranslator.hpp>
#include <glint2/HCIndex.hpp>

namespace glint2 {

void IceSheet_L0::realize()
{
	IceSheet::realize();

	// Mask out unused GCM and ice grid cells
//	HeightClassifier hc = HeightClassifier(&gcm->hcmax);
	overlap_m_hc = height_classify(
		giss::BlitzSparseMatrix(*overlap_m), elev2, gcm->hcmax);
}

/** Height Points to ice.
Uses: elev2, gcm->hpdefs, overlap */
std::unique_ptr<giss::VectorSparseMatrix> IceSheet_L0::compute_hp_to_ice()
{
	// Interpolate (for now) in height points but not X/Y
	return hp_interp(
		giss::BlitzSparseMatrix(*overlap_m), elev2, gcm->hpdefs);
}

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


/**
@param area1_m IN/OUT: Area of each GCM cell covered by
	(non-masked-out) ice sheet.
@param area1_m_hc IN/OUT: Area of each GCM cell / height class coverd by
	(non-masked-out) ice sheet
    NOTE: Indexed in 1-D according to HCIndex convention [nhc * n1] */
void IceSheet_L0::accum_areas(
giss::SparseAccumulator<int,double> &area1_m,
giss::SparseAccumulator<int,double> &area1_m_hc)
{
	accum_per_row(*overlap_m_hc, area1_m_hc);
	accum_per_row(*overlap_m, area1_m);
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
