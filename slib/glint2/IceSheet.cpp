#include <glint2/MatrixMaker.hpp>
#include <glint2/IceSheet.hpp>
#include <glint2/eigen.hpp>

namespace glint2 {

void IceSheet::clear()
{
	grid2.reset();
	exgrid.reset();
	mask2.reset();
	// elev2.clear();	// Don't know how to do this!
	overlap_raw.reset();
	overlap_m.reset();
}

/** Call this after you've set data members, to finish construction. */
void IceSheet::realize()
{
	// Compute overlap matrix from exchange grid
	overlap_raw.reset(new giss::VectorSparseMatrix(
		giss::SparseDescr(gcm->grid1->ncells_full, grid2->ncells_full)));
	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell)
		overlap_raw->add(cell->i, cell->j, area_of_polygon(*cell));

	// Mask out unused cells
	overlap_m = mask_out(
		giss::BlitzSparseMatrix(*overlap_raw), gcm->mask1.get(), mask2.get());
}

std::unique_ptr<giss::VectorSparseMatrix> IceSheet::hp_to_hc()
{
	// Get two matrices and convert to Eigen format.
	auto e_ice_to_hc(giss_to_Eigen(*ice_to_hc()));
	auto e_hp_to_ice(giss_to_Eigen(*hp_to_ice()));	// TODO: Consider making this one column-major to ease multiplication below.

	// Multiply the matices in Eigen format
	auto e_ret((*e_ice_to_hc) * (*e_hp_to_ice));

	// Convert back to GISS format sparse matrices
	return giss::Eigen_to_giss(e_ret);
}
// ==============================================================
// Write out the parts that this class computed --- so we can test/check them

boost::function<void ()> IceSheet::netcdf_define(NcFile &nc, std::string const &vname) const
{
	auto one_dim = giss::get_or_add_dim(nc, "one", 1);
	NcVar *info_var = nc.add_var((vname + ".info").c_str(), ncInt, one_dim);
	info_var->add_att("name", name().c_str());

	std::vector<boost::function<void ()>> fns;

	fns.push_back(grid2->netcdf_define(nc, vname + ".grid2"));
	fns.push_back(exgrid->netcdf_define(nc, vname + ".exgrid"));

	NcDim *n2_dim = nc.add_dim((vname + ".n2").c_str(), elev2.extent(0));
	if (mask2.get())
		fns.push_back(giss::netcdf_define(nc, vname + ".mask2", *mask2, {n2_dim}));
	fns.push_back(giss::netcdf_define(nc, vname + ".elev2", elev2, {n2_dim}));

	return boost::bind(&giss::netcdf_write_functions, fns);
}
// -------------------------------------------------------------
void IceSheet::read_from_netcdf(NcFile &nc, std::string const &vname)
{
	clear();

	grid2.reset(read_grid(nc, "grid2").release());
	exgrid.reset(read_grid(nc, "exgrid").release());
	if (giss::get_var_safe(nc, vname + ".mask2")) {
		mask2.reset(new blitz::Array<int,1>(
		giss::read_blitz<int,1>(nc, vname + ".mask2")));
	}

	elev2 = giss::read_blitz<double,1>(nc, vname + ".elev2");
}


}	// namespace glint2
