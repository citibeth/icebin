#include <glint2/MatrixMaker.hpp>

namespace glint2 {

/** See Surveyor's Formula: http://www.maa.org/pubs/Calc_articles/ma063.pdf */
inline double area_of_polygon(Cell const &cell)
{
	double ret = 0;
	auto it0 = cell.begin();
	auto it1 = it0 + 1;
	for(; it1 != cell.end(); it0=it1, it1 += 1) {
		ret += (it0->x * it1->y) - (it1->x * it0->y);
	}
	it1 = cell.begin();
	ret += (it0->x * it1->y) - (it1->x * it0->y);
	ret *= .5;
	return ret;
}

/** Computes the overlap matrix from the exchange grid */
MatrixMaker::MatrixMaker(MatrixMakerData &&data) ::
	MatrixMakerData(std::move(data))
{
	overlap_raw.reset(new giss::VectorSparseMatrix());

	for (auto cell = exgrid->cells_begin(); cell != exgrid->cells_end(); ++cell)
		overlap_raw->add(cell.i, cell.j, area_of_polygon(cell));

	overlap_m.reset(mask_out(*overlap, mask1.get(), mask2.get()));
}

std::unique_ptr<giss::VectorSparseMatrix> MatrixMaker::hp_to_hc();
{
	// Get two matrices and convert to Eigen format.
	auto e_ice_to_hc(giss_to_Eigen(*ice_to_hc()));
	auto e_hp_to_ice(giss_to_Eigen(*hp_to_ice()));	// TODO: Consider making this one column-major to ease multiplication below.

	// Multiply the matices in Eigen format
	auto e_ret((*e_ice_to_hc) * (*e_hp_to_ice));

	// Convert back to GISS format sparse matrices
	return Eigen_to_giss(*e_ret);
}



}
