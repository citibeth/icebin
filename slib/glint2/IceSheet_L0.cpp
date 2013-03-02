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
std::unique_ptr<giss::VectorSparseMatrix> IceSheet_L0::hp_to_ice()
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


// /** Computes the FHC value for grid cells relevant to this ice model.
// Does not touch FHC for other grid cells.
// @param fhc1(nhc,n1) OUT: Fraction of landice in each cell that is in the given height class.  See ModelE. */
// void IceSheet_L0::compute_fhc(
// blitz::Array<double,2> *fhc1h,
// blitz::Array<double,1> *fgice1)	// Portion of gridcell covered in ground ice (from landmask)
// {
// 	std::vector<double> area1 = overlap_raw->sum_per_row();
// 	std::vector<double> area1_m = overlap_m->sum_per_row();
// 	std::vector<double> area1_m_hc = overlap_m_hc->sum_per_row();
// 
// 	// Allow us to only iterate over cells in grid1 that overlap grid2
// 	int n1 = this->n1();
// 	giss::IndexTranslator trans_1_1p;
// 	giss::IndexTranslator trans_2_2p;
// 	make_used_translators(*overlap_m, trans_1_1p, trans_2_2p);
// 	int n1p = trans_1_1p.nb();
// 	HCIndex hci(n1);
// 	int n1hp = gcm->nhc() * n1p;
// 
// 	// Zero out ALL height classes for f1hc on any GCM grid cell that's used
// 	if (fhc1h) {
// 		// Compute f1hc
// 		// Only set items covered by this ice model
// 
// //printf("%d = (%d, %d) --> %f / %f\n", i1h, i1, hc, area1hp[i1hp], area1m[i1]);
// 		for (int i1p=0; i1p < n1p; ++i1p) {
// 			int i1 = trans_1_1p.b2a(i1p);
// 			for (int hc=0; hc<gcm->nhc(); ++hc) {
// 				int i1hc = hci.ik_to_index(i1, hc);
// 				(*fhc1h)(hc, i1) += area1_m_hc[i1hc] / area1_m[i1];
// 			}
// 		}
// 	}
// 
// 	// Only set items covered by this ice model
// 	if (fgice1)
// 	for (int i1=0; i1<n1; ++i1) {
// 		for (int i1p=0; i1p < n1p; ++i1p) {
// 			int i1 = trans_1_1p.b2a(i1p);
// 			(*fgice1)(i1) += area1_m[i1] / area1[i1];
// 		}
// 	}
// }


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



/** Computes the FHC value for grid cells relevant to this ice model.
Does not touch FHC for other grid cells.
@param fhc1(nhc,n1) OUT: Fraction of landice in each cell that is in the given height class.  See ModelE. */
void IceSheet_L0::compute_fhc2(
giss::SparseAccumulator<int,double> &fgice1,
giss::SparseAccumulator<int,double> &fhc1h,
giss::SparseAccumulator<int,double> &area1_m)
{
	// ======== Compute fhc1h (must divide by area1_m afterwards)
	SparseAccumulator<int,double> area1_m_hc;
	accum_per_row(*overlap_m_hc, area1_m_hc);



	HCIndex hci(this->n1());	// Calculate bundled hc/i1 indices



		for (int hc=0; hc<gcm->nhc(); ++hc) {
			int i1hc = hci.ik_to_index(i1, hc);
			auto ii(area1_m_hc.find(i1hc));
			if (ii != area1_m_hc.end()) {
				// fhc1h_vals.push_back(ii->second / area1_m[i1]);
				fhc1h_vals.push_back(ii->second);
			} else fhc1h_vals.push_back(0);




	SparseAccumulator<int,double> area1_m;
	accum_per_row(*overlap_m, area1_m);



	// Compute fgice1
	for (auto ii = area1_m.begin(); ii != area1_m.end(); ++i) {
		int i1 = ii->first;
		double area1_m_i1 = ii->second;

		indices.push_back(i1);
		fgice1_vals.push_back(area1_m_i1 / gcm->grid1->get_cell(i1)->area);
	}

	// Compute fhc1h
	HCIndex hci(this->n1());	// Calculate bundled hc/i1 indices



		for (int hc=0; hc<gcm->nhc(); ++hc) {
			int i1hc = hci.ik_to_index(i1, hc);
			auto ii(area1_m_hc.find(i1hc));
			if (ii != area1_m_hc.end()) {
				// fhc1h_vals.push_back(ii->second / area1_m[i1]);
				fhc1h_vals.push_back(ii->second);
			} else fhc1h_vals.push_back(0);
	}



	std::vector<double> area1_m = overlap_m->sum_per_row();
	std::vector<double> area1_m_hc = overlap_m_hc->sum_per_row();

	// Allow us to only iterate over cells in grid1 that overlap grid2
	int n1 = this->n1();
	giss::IndexTranslator trans_1_1p;
	giss::IndexTranslator trans_2_2p;
	make_used_translators(*overlap_m, trans_1_1p, trans_2_2p);
	int n1p = trans_1_1p.nb();
	HCIndex hci(n1);
	int n1hp = gcm->nhc() * n1p;

	// Compute f1hc
	// Only set items covered by this ice model

//printf("%d = (%d, %d) --> %f / %f\n", i1h, i1, hc, area1hp[i1hp], area1m[i1]);
	for (int i1p=0; i1p < n1p; ++i1p) {
		int i1 = trans_1_1p.b2a(i1p);

		indices1.push_back(i1);

		// Me shoudl be dividing by area of the gridcell, not just the overlap area
		// fgice1_vals.push_back(area1_m[i1] / area1[i1]);
		fgice1_vals.push_back(area1_m[i1] /
			gcm->grid1->get_cell(i1)->area);

		for (int hc=0; hc<gcm->nhc(); ++hc) {
			int i1hc = hci.ik_to_index(i1, hc);
			fhc1h_vals.push_back(area1_m_hc[i1hc] / area1_m[i1]);
...for multiple ice sheets, we need to be dividing by sum(area1_m[i1]) over all ice sheets...
			// (*fhc1h)(hc, i1) += area1_m_hc[i1hc] / area1_m[i1];
		}
	}
}




boost::function<void ()> IceSheet_L0::netcdf_define(NcFile &nc, std::string const &vname) const
{
	auto ret = IceSheet::netcdf_define(nc, vname);

	NcVar *info_var = nc.get_var((vname + ".info").c_str());
	info_var->add_att("parameterization", "L0");

	return ret;
}
// -------------------------------------------------------------


}	// namespace glint2
