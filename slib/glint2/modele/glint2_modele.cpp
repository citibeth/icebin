#include <netcdfcpp.h>
#include <giss/blitz.hpp>
#include <giss/f90blitz.hpp>
#include <glint2/HCIndex.hpp>
#include <glint2/modele/glint2_modele.hpp>

using namespace glint2;
using namespace glint2::modele;


/** @param spec  */
extern "C" glint2_modele *glint2_modele_new(
	char const *maker_fname_f, int maker_fname_len,
	char const *maker_vname_f, int maker_vname_len,

	// Info about the global grid
	int im, int jm,

	// Info about the local grid (C-style indices)
	int i0h, int i1h, int j0h, int j1h,
	int i0, int i1, int j0, int j1,
	int j0s, int j1s,

	// MPI Stuff
	int comm_f, int root)
{
	printf("***** BEGIN glint2_modele()\n");

	// Convert Fortran arguments
	std::string maker_fname(maker_fname_f, maker_fname_len);
	std::string maker_vname(maker_vname_f, maker_vname_len);

	// Allocate our return variable
	std::unique_ptr<glint2_modele> api(new glint2_modele());

	// Set up the domain
	std::unique_ptr<GridDomain> mdomain(
		new ModelEDomain(im, jm,
			i0h, i1h, j0h, j1h,
			i0, i1, j0, j1,
			j0s, j1s));
	api->domain = (ModelEDomain *)mdomain.get();

	// Load the MatrixMaker	(filtering by our domain, of course)
	// Also load the ice sheets
	api->maker.reset(new MatrixMaker(std::move(mdomain)));
	NcFile nc(maker_fname.c_str(), NcFile::ReadOnly);
	api->maker->read_from_netcdf(nc, maker_vname);
	api->maker->realize();

	// Read the coupler, along with ice model proxies
	api->gcm_coupler.reset(new GCMCoupler_MPI(MPI_Comm_f2c(comm_f), root));
	api->gcm_coupler->read_from_netcdf(nc, maker_vname, api->maker->get_sheet_names(), api->maker->sheets);
	nc.close();



	// TODO: Test that im and jm are consistent with the grid read.

	printf("***** END glint2_modele()\n");

	// No exception, we can release our pointer back to Fortran
	return api.release();
}
// -----------------------------------------------------
extern "C" void glint2_modele_delete(glint2_modele *&api)
{
	if (api) delete api;
	api = 0;
}
// -----------------------------------------------------
extern "C"
int glint2_modele_nhc(glint2_modele *api)
{
	int ret = api->maker->nhc();
	printf("glint2_modele_nhc() returning %d\n", ret);
	return ret;
}
// -----------------------------------------------------
extern "C"
void glint2_modele_get_elevhc_c(glint2_modele *api,
	giss::F90Array<double, 3> &elevhc_f)			// OUT
{
	auto elevhc(elevhc_f.to_blitz());
	int nhc = api->maker->nhc();
	if (nhc != elevhc.extent(2)) {
		fprintf(stderr, "glint2_modele_get_elevhc: Inconsistent nhc (%d vs %d)\n", elevhc.extent(2), api->maker->nhc());
		throw std::exception();
	}

	for (int k=1; k <= nhc; ++k) {
		double val = api->maker->hpdefs[k-1];
		for (int j=elevhc.lbound(1); j <= elevhc.ubound(1); ++j) {
		for (int i=elevhc.lbound(0); i <= elevhc.ubound(0); ++i) {
			elevhc(i,j,k) = val;
		}
	}}
}
// -----------------------------------------------------
extern "C"
void glint2_modele_compute_fgice_c(glint2_modele *api,
	giss::F90Array<double, 2> &fgice1_f,		// OUT
	giss::F90Array<double, 2> &fgrnd1_f,		// IN/OUT
	giss::F90Array<double, 2> &focean1_f,		// IN
	giss::F90Array<double, 2> &flake1_f			// IN
)
{
	ModelEDomain &domain(*api->domain);

#if 0
printf("domain: (%d %d) (%d %d %d %d) (%d %d %d %d) (%d %d)\n",
domain.im, domain.jm,
domain.i0h_f, domain.i1h_f, domain.j0h_f, domain.j1h_f,
domain.i0_f, domain.i1_f, domain.j0_f, domain.j1_f,
domain.j0s_f, domain.j1s_f);

printf("grid1->size() = %ld\n", api->maker->grid1->ncells_realized());
#endif

	// Reconstruct arrays, using Fortran conventions
	// (smallest stride first, whatever-based indexing it came with)

	// Get the sparse vector values
	giss::CooVector<std::pair<int,int>,double> fhc1h_s;
	giss::CooVector<int,double> fgice1_s;
	api->maker->compute_fgice(fgice1_s);

	// Translate the sparse vectors to the ModelE data structures
	std::vector<std::tuple<int, int, double>> fgice1_vals;
	for (auto ii = fgice1_s.begin(); ii != fgice1_s.end(); ++ii) {
		int i1 = ii->first;

		// Filter out things not in our domain
		// (we'll get the answer for our halo via a halo update)
		// Convert to local (ModelE 2-D) indexing convention
		int lindex[domain.num_local_indices];
		domain.global_to_local(i1, lindex);
		if (!domain.in_domain(lindex)) continue;

		// Store it away
		// (we've eliminated duplicates, so += isn't needed, but doesn't hurt either)
		fgice1_vals.push_back(std::make_tuple(lindex[0], lindex[1], ii->second));
	}

	// Zero out fgice1, ONLY where we're touching it.
	auto fgice1(fgice1_f.to_blitz());
	for (auto ii=fgice1_vals.begin(); ii != fgice1_vals.end(); ++ii) {
		int ix_i = std::get<0>(*ii);
		int ix_j = std::get<1>(*ii);
		double val = std::get<2>(*ii);

		fgice1(ix_i, ix_j) = 0;
	}

	// Replace with our values
	for (auto ii=fgice1_vals.begin(); ii != fgice1_vals.end(); ++ii) {
		int ix_i = std::get<0>(*ii);
		int ix_j = std::get<1>(*ii);
		double val = std::get<2>(*ii);

		fgice1(ix_i, ix_j) += val;
	}
	// -----------------------------------------------------
	// Balance fgice against other landcover types
	auto fgrnd1(fgrnd1_f.to_blitz());
	auto focean1(focean1_f.to_blitz());
	auto flake1(flake1_f.to_blitz());
	fgrnd1 = 1.0 - focean1 - flake1 - fgice1;

	// -----------------------------------------------------
}
// -----------------------------------------------------
extern "C"
void glint2_modele_compute_fhc_c(glint2_modele *api,
	giss::F90Array<double, 3> &fhc1h_f)			// IN/OUT
{
	ModelEDomain &domain(*api->domain);

#if 0
printf("domain: (%d %d) (%d %d %d %d) (%d %d %d %d) (%d %d)\n",
domain.im, domain.jm,
domain.i0h_f, domain.i1h_f, domain.j0h_f, domain.j1h_f,
domain.i0_f, domain.i1_f, domain.j0_f, domain.j1_f,
domain.j0s_f, domain.j1s_f);

printf("grid1->size() = %ld\n", api->maker->grid1->ncells_realized());
#endif

	// Reconstruct arrays, using Fortran conventions
	// (smallest stride first, whatever-based indexing it came with)

	// Get the sparse vector values
	giss::CooVector<std::pair<int,int>,double> fhc1h_s;
	api->maker->compute_fhc(fhc1h_s);

	std::vector<std::tuple<int, int, int, double>> fhc1h_vals;

	// Work on fhc1h
	for (auto ii = fhc1h_s.begin(); ii != fhc1h_s.end(); ++ii) {
		int i1 = ii->first.first;
		int hc = ii->first.second;		// zero-based
		double val = ii->second;

		// Filter out things not in our domain
		// (we'll get the answer for our halo via a halo update)

		// Convert to local (ModelE 2-D + height class) indexing convention
		int lindex[domain.num_local_indices];
		domain.global_to_local(i1, lindex);
		if (!domain.in_domain(lindex)) continue;

		// Store it away
		// (we've eliminated duplicates, so += isn't needed, but doesn't hurt either)
		int hc_f = hc + 1;		// convert zero-based to 1-based arrays
		fhc1h_vals.push_back(std::make_tuple(lindex[0], lindex[1], hc_f, val));
	}

	auto fhc1h(fhc1h_f.to_blitz());
	fhc1h = 0;
	fhc1h(blitz::Range::all(), blitz::Range::all(), 1) = 1.0;
	for (auto ii=fhc1h_vals.begin(); ii != fhc1h_vals.end(); ++ii) {
		int ix_i = std::get<0>(*ii);
		int ix_j = std::get<1>(*ii);
		int hc_f = std::get<2>(*ii);
		double val = std::get<3>(*ii);

		fhc1h(ix_i, ix_j, 1) = 0;
	}

	for (auto ii=fhc1h_vals.begin(); ii != fhc1h_vals.end(); ++ii) {
		int ix_i = std::get<0>(*ii);
		int ix_j = std::get<1>(*ii);
		int hc_f = std::get<2>(*ii);
		double val = std::get<3>(*ii);

		fhc1h(ix_i, ix_j, hc_f) += val;
	}
}
// -----------------------------------------------------



static std::unique_ptr<giss::VectorSparseMatrix> filter_matrix_hp(
	HCIndex const &hc_index,
	GridDomain const &domain1,
	GridDomain const &domain2,
	giss::VectorSparseMatrix const &mat)
{
	std::unique_ptr<giss::VectorSparseMatrix> ret(
		new giss::VectorSparseMatrix((giss::SparseDescr)mat));

	int lindex1[domain1.num_local_indices];
	int lindex2[domain2.num_local_indices];
	for (auto ii = mat.begin(); ii != mat.end(); ++ii) {

		// Output of linear transformation: Only include
		// if it's part of our domain.
		int hc1, i1;
		hc_index.index_to_ik(ii.row(), i1, hc1);
		domain1.global_to_local(i1, lindex1);
		if (!domain1.in_domain(lindex1)) {
//printf("Throwing out of domain: %d %d\n", lindex1[0], lindex1[1]);
			continue;
		}

		// Input of linear transformation: must be in halo
		int hc2, i2;
		hc_index.index_to_ik(ii.col(), i2, hc2);
		domain2.global_to_local(i2, lindex2);
		if (!domain2.in_halo(lindex2)) {
			fprintf(stderr, "Error filtering matrix: grid cell %d (", ii.col());
			for (int i=0; i<domain2.num_local_indices; ++i) fprintf(stderr, "%d ", lindex2[i]);
			fprintf(stderr, ") in input (column) is not available in the halo.\n");
			throw std::exception();
		}

		ret->add(ii.row(), ii.col(), ii.val());
	}

printf("filter_matrix_hp went from size %ld to %ld\n", mat.size(), ret->size());
	return ret;
}



/** Call this to figure out how to dimension arrays.
@return Number of elements in the sparse matrix */
extern "C"
int glint2_modele_hp_to_hc_part1(glint2_modele *api)
{
	auto mat(api->maker->hp_to_hc());
	HCIndex hc_index(api->maker->n1());
	api->hp_to_hc = filter_matrix_hp(hc_index, *api->domain, *api->domain, *mat);
printf("Filtered matrix at %p (from %p)\n", api->hp_to_hc.get(), mat.get());
	// api->hp_to_hc = std::move(mat);	// debugging
	return api->hp_to_hc->size();
}
// -----------------------------------------------------
static void global_to_local_hp(
	glint2_modele *api,	
	HCIndex const &hc_index,
	std::vector<int> const &grows,
	blitz::Array<int,1> rows_i,		// Fortran-style array, base=1
	blitz::Array<int,1> rows_j,
	blitz::Array<int,1> rows_k)		// height point index
{
	// Copy the rows while translating
	// auto rows_k(rows_k_f.to_blitz());
	//std::vector<double> &grows = *api->hp_to_hc.rows();
	int lindex[api->domain->num_local_indices];
	for (int i=0; i<grows.size(); ++i) {		
		int ihc, i1;
		hc_index.index_to_ik(grows[i], i1, ihc);
		api->domain->global_to_local(i1, lindex);
//if (lindex[1] >= 80) printf("Found big lindex: %d %d\n", lindex[0], lindex[1]);
		rows_i(i+1) = lindex[0];
		rows_j(i+1) = lindex[1];
		rows_k(i+1) = ihc+1;	// Convert to Fortran indexing
	}
}
// -----------------------------------------------------
/** Call this after rows, cols and vals have been dimensioned. */
extern "C"
void glint2_modele_hp_to_hc_part2(glint2_modele *api,
	giss::F90Array<int, 1> &rows_i_f,
	giss::F90Array<int, 1> &rows_j_f,
	giss::F90Array<int, 1> &rows_k_f,
	giss::F90Array<int, 1> &cols_i_f,
	giss::F90Array<int, 1> &cols_j_f,
	giss::F90Array<int, 1> &cols_k_f,
	giss::F90Array<double, 1> &vals_f)
{

	// Array bounds checking not needed, it's done
	// in lower-level subroutines that we call.

	HCIndex hc_index(api->maker->n1());

printf("Translating matrix at %p\n", api->hp_to_hc.get());
	// Translate rows and cols
	global_to_local_hp(api, hc_index, api->hp_to_hc->rows(),
		rows_i_f.to_blitz(),
		rows_j_f.to_blitz(),
		rows_k_f.to_blitz());
printf("Translating: Done With Rows!\n");
	global_to_local_hp(api, hc_index, api->hp_to_hc->cols(),
		cols_i_f.to_blitz(),
		cols_j_f.to_blitz(),
		cols_k_f.to_blitz());

	// Copy the values, just a simple vector copy
	auto vals(vals_f.to_blitz());
	std::vector<double> const &mvals(api->hp_to_hc->vals());
	for (int i=0; i<mvals.size(); ++i) vals(i+1) = mvals[i];

	// Free temporary storage
	api->hp_to_hc.reset();
}
// -----------------------------------------------------
/** @param hpvals Values on height-points GCM grid for various fields
	the GCM has decided to provide. */
extern "C"
void glint2_modele_couple_to_ice(
glint2_modele *api,
giss::F90Array<double,3> &smb1hp_f,
giss::F90Array<double,3> &seb1hp_f)
{
	std::vector<IceField> fields =
		{IceField::MASS_FLUX, IceField::ENERGY_FLUX};
//	std::vector<blitz::Array<double,3>> vals1hp =
//		{smb1hp_f.to_blitz(), seb1hp_f.to_blitz()};

	auto smb1hp(smb1hp_f.to_blitz());
	auto seb1hp(seb1hp_f.to_blitz());

	// Count total number of elements in the matrices
	// (_l = local to this MPI node)
	int nele_l = api->maker->ice_matrices_size();

	// Allocate buffer for that amount of stuff
	int nfields = fields.size();
	giss::DynArray<SMBMsg> sbuf(SMBMsg::size(nfields), nele_l);

	// Fill it in by doing a sparse multiply...
	// (while translating indices to local coordinates)
	HCIndex hc_index(api->maker->n1());
	int nmsg = 0;
	for (auto sheet=api->maker->sheets.begin(); sheet != api->maker->sheets.end(); ++sheet) {
		int sheetno = sheet->index;
//printf("glint2_modele: %p.sheetno = %d\n", &*sheet, sheetno);
		giss::VectorSparseMatrix &mat(sheet->hp_to_ice());

		// Skip if we have nothing to do for this ice sheet
		if (mat.size() == 0) continue;

		// Do the multiplication
		for (auto ii=mat.begin(); ii != mat.end(); ++ii) {
			SMBMsg &msg = sbuf[nmsg];
			msg.sheetno = sheetno;
			msg.i2 = ii.row();

			int i1, ihc;
			hc_index.index_to_ik(ii.col(), i1, ihc);
			int lindex[api->domain->num_local_indices];
			api->domain->global_to_local(i1, lindex);
			msg[0] = ii.val() * smb1hp(lindex[0], lindex[1], ihc+1);
			msg[1] = ii.val() * seb1hp(lindex[0], lindex[1], ihc+1);
//printf("msg = %d (i,j, hc)=(%d %d %d) i2=%d %g %g (%g %g)\n", msg.sheetno, lindex[0], lindex[1], ihc+1, msg.i2, msg[0], msg[1], smb1hp(lindex[0], lindex[1], ihc+1), seb1hp(lindex[0], lindex[1], ihc+1));

			++nmsg;
		}
	}

	// Sanity check: make sure we haven't overrun our buffer
	if (nmsg != sbuf.size) {
		fprintf(stderr, "Wrong number of items in buffer: %d vs %d expected\n", nmsg, sbuf.size);
		throw std::exception();
	}

	api->gcm_coupler->couple_to_ice(fields, sbuf);
}
