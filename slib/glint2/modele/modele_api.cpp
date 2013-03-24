#include <netcdfcpp.h>
#include <giss/blitz.hpp>
#include <giss/f90blitz.hpp>
#include <glint2/HCIndex.hpp>
#include <glint2/modele_api.hpp>

using namespace glint2;
using namespace glint2::modele;


/** @param spec  */
extern "C" modele_api *modele_api_new(
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
	// Convert Fortran arguments
	std::string maker_fname(maker_fname_f, maker_fname_len);
	std::string maker_vname(maker_vname_f, maker_vname_len);

	// Allocate our return variable
	std::unique_ptr<modele_api> api(new modele_api());

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

	// No exception, we can release our pointer back to Fortran
	return api.release();
}
// -----------------------------------------------------
extern "C" void modele_api_delete(modele_api *&api)
{
	if (api) delete api;
	api = 0;
}
// -----------------------------------------------------
extern "C"
void modele_api_compute_fhc_c(modele_api *api,
	giss::F90Array<double, 3> &fhc1h_f,			// IN/OUT
	giss::F90Array<double, 2> &fgice1_f,		// IN/OUT
	giss::F90Array<double, 2> &fgrnd1_f,		// OUT
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
	api->maker->compute_fhc(fhc1h_s, fgice1_s);

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
/** Call this to figure out how to dimension arrays.
@return Number of elements in the sparse matrix */
extern "C"
int modele_api_hp_to_hc_part1(modele_api *api)
{
	auto mat(api->maker->hp_to_hc());
	api->hp_to_hc = filter_matrix(*api->domain, *api->domain, *mat);
	return api->hp_to_hc->size();
}
// -----------------------------------------------------
static void global_to_local_hp(
	modele_api *api,	
	HCIndex const &hc_index,
	std::vector<int> const &grows,
	blitz::Array<int,1> rows_i,
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
		rows_i(i) = lindex[0];
		rows_j(i) = lindex[1];
		rows_k(i) = ihc+1;	// Convert to Fortran indexing
	}
}
// -----------------------------------------------------
/** Call this after rows, cols and vals have been dimensioned. */
extern "C"
void modele_api_hp_to_hc_part2(modele_api *api,
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

	// Translate rows and cols
	global_to_local_hp(api, hc_index, api->hp_to_hc->rows(),
		rows_i_f.to_blitz(),
		rows_j_f.to_blitz(),
		rows_k_f.to_blitz());
	global_to_local_hp(api, hc_index, api->hp_to_hc->cols(),
		cols_i_f.to_blitz(),
		cols_j_f.to_blitz(),
		cols_k_f.to_blitz());

	// Copy the values, just a simple vector copy
	auto vals(vals_f.to_blitz());
	vals = giss::vector_to_blitz(api->hp_to_hc->vals());

	// Free temporary storage
	api->hp_to_hc.reset();
}
// -----------------------------------------------------
/** @param hpvals Values on height-points GCM grid for various fields
	the GCM has decided to provide. */
extern "C"
void modele_api_couple_to_ice(
modele_api *api,
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
//printf("modele_api: %p.sheetno = %d\n", &*sheet, sheetno);
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
