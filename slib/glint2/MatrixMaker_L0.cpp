#include <glint2/MatrixMaker.hpp>

namespace glint2 {

MatrixMaker_L0::MatrixMaker_L0(MatrixMakerData &&data) : MatrixMaker(std::move(data))
{
	// Mask out unused GCM and ice grid cells
	overmap_m_hc.reset(
		height_classify(*overlap_m, elev2, HeightClassifier(&hcmax)));
}

/** Height Points to ice.
Uses: elev2, hpdefs, overlap */
std::unique_ptr<VectorSparseMatrix> MatrixMaker_L0::hp_to_ice()
{
	// Interpolate (for now) in height points but not X/Y
	return hc_interp(*overlap_m, hpdefs, elev2);
}

/** Uses: elev2, hcmax, overlap */
virtual std::unique_ptr<VectorSparseMatrix> ice_to_hc()
{
	// Just area-weighted remap from ice back to hc grid
	auto g2_to_g1(grid2_to_grid1(overlap_m_hc));
	proj_to_native(grid1, *g2_to_g1);

	return g2_to_g1;
}



/** Computes the FHC value for grid cells relevant to this ice model.
Does not touch FHC for other grid cells.
@param fhc1(nhc,n1) OUT: Fraction of landice in each cell that is in the given height class.  See ModelE. */
void MatrixMaker_L0::compute_fhc(
blitz::Array<double,2> *fhc1h,
blitz::Array<double,2> *elev1h,
blitz::Array<double,1> *fgice1)	// Portion of gridcell covered in ground ice (from landmask)
{
	std::vector<double> area1 = overlap_raw->sum_per_row();
	std::vector<double> area1_m = overlap_m->sum_per_row();
	std::vector<double> area1_m_hc = overlap_m_hc->sum_per_row();

	// Allow us to only iterate over cells in grid1 that overlap grid2
	IndexTranslator trans_1_1p;
	IndexTranslator trans_2_2p;
	make_used_translators(*overlap_m, trans_1_1p, trans_2_2p);
	int n1p = trans_1_1p.nb();

	// Zero out ALL height classes for f1hc on any GCM grid cell that's used
	if (fhc1h) {
		std::vector<bool> zerome(n1,false);
		for (int i1p=0; i1p < n1p; ++i1p) {
			int i1 = trans_1_1p.b2a(i1p);

			for (int ihc=0; ihc<nhc; ++ihc)
				(*fhc1h)(ihc, i1) = 0;
		}
	}

	// Compute f1hc
	// Only set items covered by this ice model
	for (int i1hp=0; i1hp < n1hp; ++i1hp) {
		int i1h = i1hp_to_i1h(i1hp);
		int i1 = i1h_to_i1(i1h);
		int hc = get_hclass(i1h, i1);

//printf("%d = (%d, %d) --> %f / %f\n", i1h, i1, hc, area1hp[i1hp], area1m[i1]);
		if (fhc1h) (*fhc1h)(hc, i1) = area1hp[i1hp] / area1m[i1];
		if (elev1h) (*elev1h)(hc, i1) = elev1hp[i1hp];
	}

	// Only set items covered by this ice model
	if (fgice1)
	for (int i1=0; i1<n1; ++i1) {
		if (area1m[i1] > 0)
			(*fgice1)(i1) = area1m[i1] / area1[i1];
	}
}

}	// namespace glint2
