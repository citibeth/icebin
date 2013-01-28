#include <glint2/MatrixMaker.hpp>
#include <glint2/eigen.hpp>

namespace glint2 {

/** Call this after you've set data members, to finish construction. */
void MatrixMaker::realize()
{

	// Compute overlap matrix from exchange grid
	overlap_raw.reset(new giss::VectorSparseMatrix(
		giss::SparseDescr(grid1->ncells_full, grid2->ncells_full)));
	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell)
		overlap_raw->add(cell->i, cell->j, area_of_polygon(*cell));

	// Mask out unused cells
	overlap_m = mask_out(
		giss::BlitzSparseMatrix(*overlap_raw), mask1.get(), mask2.get());
}

std::unique_ptr<giss::VectorSparseMatrix> MatrixMaker::hp_to_hc()
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

static void MatrixMaker_netcdf_write(// NcFile *nc, std::string const &vname,
boost::function<void ()> const &fn1,
boost::function<void ()> const &fn2)
{
	fn1();
	fn2();
}


boost::function<void ()> MatrixMaker::netcdf_define(NcFile &nc, std::string const &vname) const
{
	// ------ Attributes
	auto one_dim = giss::get_or_add_dim(nc, "one", 1);
	NcVar *info_var = nc.add_var((vname + ".info").c_str(), ncInt, one_dim);
		info_var->add_att("grid1.name", grid1->name.c_str());
		info_var->add_att("grid2.name", grid2->name.c_str());
		info_var->add_att("exgrid.name", exgrid->name.c_str());

	return boost::bind(&MatrixMaker_netcdf_write, // this, &nc, vname,
		overlap_raw->netcdf_define(nc, vname + ".overlap_raw"),
		overlap_m->netcdf_define(nc, vname + ".overlap_m"));
}




}
