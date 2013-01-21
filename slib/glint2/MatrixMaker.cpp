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
	overlap_m = mask_out(*overlap_raw, mask1.get(), mask2.get());
}

std::unique_ptr<giss::VectorSparseMatrix> MatrixMaker::hp_to_hc()
{
	// Get two matrices and convert to Eigen format.
	auto e_ice_to_hc(giss_to_Eigen(*ice_to_hc()));
	auto e_hp_to_ice(giss_to_Eigen(*hp_to_ice()));	// TODO: Consider making this one column-major to ease multiplication below.

	// Multiply the matices in Eigen format
	auto e_ret((*e_ice_to_hc) * (*e_hp_to_ice));

	// Convert back to GISS format sparse matrices
	return Eigen_to_giss(e_ret);
}



}
