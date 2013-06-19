#include <mpi.h>	// For Intel MPI, mpi.h must be included before stdio.h
#include <netcdfcpp.h>
#include <giss/blitz.hpp>
#include <giss/f90blitz.hpp>
#include <glint2/HCIndex.hpp>
#include <glint2/modele/glint2_modele.hpp>
#include <glint2/IceModel_TConv.hpp>

using namespace glint2;
using namespace glint2::modele;

// ---------------------------------------------------

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
	int comm_f, int root,

	// Constants from ModelE
	double LHM, double SHI)
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

	// Add adapters to the ice models, if needed
	for (auto ii=api->gcm_coupler->models.super::begin();
		ii != api->gcm_coupler->models.super::end(); ++ii)
	{
		std::set<IceField> fields;
		ii->second->get_required_fields(fields);
		if (fields.find(IceField::SURFACE_T) != fields.end()) {
			IceModel_Decode *imd = dynamic_cast<IceModel_Decode *>(&*(ii->second));
			if (!imd) continue;

			// This ice model wants SURFACE_T... add an adapter
			ii->second.release();
			std::unique_ptr<IceModel> model(std::move(ii->second));
			ii->second.reset(new IceModel_TConv(
				std::unique_ptr<IceModel_Decode>(imd),
				LHM, SHI));
		}
	}


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
int glint2_modele_nhp(glint2_modele *api)
{
	int ret = api->maker->nhp();
	// HP/HC = 1 (Fortran) reserved for legacy "non-model" ice
    // (not part of GLINT2)
	ret += 1;
	printf("glint2_modele_nhp() returning %d\n", ret);
	return ret;
}
// -----------------------------------------------------
extern "C"
void glint2_modele_compute_fgice_c(glint2_modele *api,
	int replace_fgice_b,
	giss::F90Array<double, 2> &fgice1_glint2_f,		// OUT
	giss::F90Array<double, 2> &fgice1_f,		// OUT
	giss::F90Array<double, 2> &fgrnd1_f,		// IN/OUT
	giss::F90Array<double, 2> &focean1_f,		// IN
	giss::F90Array<double, 2> &flake1_f			// IN
)
{
printf("BEGIN glint2_modele_compute_fgice_c()\n");
	ModelEDomain &domain(*api->domain);

#if 0
{
NcFile nc("fgice1_0.nc", NcFile::Replace);
auto fgice1(fgice1_f.to_blitz());
auto fgice1_c(giss::f_to_c(fgice1));
auto a(giss::netcdf_define(nc, "fgice1", fgice1_c));
a();
nc.close();
}
#endif

#if 0
printf("domain: (%d %d) (%d %d %d %d) (%d %d %d %d) (%d %d)\n",
domain.im, domain.jm,
domain.i0h_f, domain.i1h_f, domain.j0h_f, domain.j1h_f,
domain.i0_f, domain.i1_f, domain.j0_f, domain.j1_f,
domain.j0s_f, domain.j1s_f);

printf("grid1->size() = %ld\n", api->maker->grid1->ncells_realized());
#endif

//printf("Addresses:\n%p\n%p\n%p\n%p\n%p\n", fgice1_glint2_f.base, fgice1_f.base, fgrnd1_f.base, focean1_f.base, flake1_f.base);

	// Reconstruct arrays, using Fortran conventions
	// (smallest stride first, whatever-based indexing it came with)

	// Get the sparse vector values
	giss::CooVector<std::pair<int,int>,double> fhc1h_s;
	giss::CooVector<int,double> fgice1_s;
	api->maker->fgice(fgice1_s);

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
	if (replace_fgice_b != 0) {
		for (auto ii=fgice1_vals.begin(); ii != fgice1_vals.end(); ++ii) {
			int ix_i = std::get<0>(*ii);
			int ix_j = std::get<1>(*ii);
			// double val = std::get<2>(*ii);

			fgice1(ix_i, ix_j) = 0;
		}
	}

	// Zero out the GLINT2-only version completely
	auto fgice1_glint2(fgice1_glint2_f.to_blitz());
	fgice1_glint2 = 0;

	// Replace with our values
	for (auto ii=fgice1_vals.begin(); ii != fgice1_vals.end(); ++ii) {
		int ix_i = std::get<0>(*ii);
		int ix_j = std::get<1>(*ii);
		double val = std::get<2>(*ii);

		fgice1(ix_i, ix_j) += val;
		fgice1_glint2(ix_i, ix_j) += val;
	}
	// -----------------------------------------------------
	// Balance fgice against other landcover types
	auto fgrnd1(fgrnd1_f.to_blitz());
	auto focean1(focean1_f.to_blitz());
	auto flake1(flake1_f.to_blitz());
	fgrnd1 = 1.0 - focean1 - flake1 - fgice1;

#if 0
printf("AA1\n");
printf("fgrnd1: (%d %d) - (%d %d) [%d %d]\n",
	fgrnd1.lbound(0), fgrnd1.lbound(1),
	fgrnd1.ubound(0), fgrnd1.ubound(1),
	fgrnd1.extent(0), fgrnd1.extent(1));

NcFile nc("fgice1_1.nc", NcFile::Replace);
auto fgice1_c(giss::f_to_c(fgice1));
auto fgice1_glint2_c(giss::f_to_c(fgice1_glint2));
auto a(giss::netcdf_define(nc, "fgice1", fgice1_c));
auto b(giss::netcdf_define(nc, "fgice1_glint2", fgice1_glint2_c));
a();
b();
nc.close();
#endif

	// -----------------------------------------------------
printf("END glint2_modele_compute_fgice_c()\n");
}
// -----------------------------------------------------
static void global_to_local_hp(
	glint2_modele *api,	
	HCIndex const &hc_index,
	std::vector<int> const &grows,
	std::string const &name,	// For debugging
	blitz::Array<int,1> &rows_i,		// Fortran-style array, base=1
	blitz::Array<int,1> &rows_j,
	blitz::Array<int,1> &rows_k)		// height point index
{
printf("BEGIN global_to_local_hp %p %p %p %p\n", &grows[0], rows_i.data(), rows_j.data(), rows_k.data());
	// Copy the rows while translating
	// auto rows_k(rows_k_f.to_blitz());
	//std::vector<double> &grows = *api->hp_to_hc.rows();
	int lindex[api->domain->num_local_indices];
	for (int i=0; i<grows.size(); ++i) {		
		int ihc, i1;
		hc_index.index_to_ik(grows[i], i1, ihc);
		api->domain->global_to_local(i1, lindex);
		rows_i(i+1) = lindex[0];
		rows_j(i+1) = lindex[1];
		// +1 for C-to-Fortran conversion
		// +1 because lowest HP/HC is reserved
		rows_k(i+1) = ihc+2;
	}
printf("END global_to_local_hp\n");
}
// -----------------------------------------------------
/**
@param zatmo1_f ZATMO from ModelE (Elevation of bottom of atmosphere * GRAV)
@param BYGRAV 1/GRAV = 1/(9.8 m/s^2)
@param fgice1_glint2_f Amount of GLINT2-related ground ice in each GCM grid cell
@param fgice1_f Total amount of ground ice in each GCM grid cell
@param used1h_f Height point mask
@param fhc1h_f Weights to average height points into GCM grid.
@param elev1h_f Elevation of each height point
*/
extern "C"
void glint2_modele_init_landice_com_c(glint2::modele::glint2_modele *api,
	giss::F90Array<double, 2> &zatmo1_f,	// IN
	double const BYGRAV,					// IN
	giss::F90Array<double, 2> &fgice1_glint2_f,	// IN
	giss::F90Array<double, 2> &fgice1_f,	// IN
	giss::F90Array<int,3> &used1h_f,		// IN/OUT
	giss::F90Array<double, 3> &fhc1h_f,		// OUT: hp-to-atmosphere
	giss::F90Array<double, 3> &elev1h_f)	// IN/OUT
{
printf("init_landice_com_part2 1\n");

	// =================== elev1h
	// Just copy out of hpdefs array, elevation points are the same
	// on all grid cells.

	auto elev1h(elev1h_f.to_blitz());
	int nhp_glint2 = api->maker->nhp();
	int nhp = api->maker->nhp() + 1;	// Add non-model HP
	if (nhp != elev1h.extent(2)) {
		fprintf(stderr, "glint2_modele_get_elev1h: Inconsistent nhp (%d vs %d)\n", elev1h.extent(2), nhp);
		throw std::exception();
	}

	// Copy 1-D height point definitions to elev1h
	for (int k=0; k < nhp_glint2; ++k) {
		double val = api->maker->hpdefs[k];
		for (int j=elev1h.lbound(1); j <= elev1h.ubound(1); ++j) {
		for (int i=elev1h.lbound(0); i <= elev1h.ubound(0); ++i) {
			// +1 for C-to-Fortran conversion
			// +1 because lowest HP/HC is reserved
			elev1h(i,j,k+2) = val;
		}}
	}

	// Copy zatmo to elevation of reserved height point
	auto zatmo1(zatmo1_f.to_blitz());
	for (int j=elev1h.lbound(1); j <= elev1h.ubound(1); ++j) {
	for (int i=elev1h.lbound(0); i <= elev1h.ubound(0); ++i) {
		elev1h(i,j,1) = zatmo1(i,j) * BYGRAV;
	}}

printf("init_landice_com_part2 2\n");
	// ======================= fhc(:,:,1)
	auto fgice1(fgice1_f.to_blitz());
	auto fgice1_glint2(fgice1_glint2_f.to_blitz());
	auto fhc1h(fhc1h_f.to_blitz());
	fhc1h = 0;
	for (int j=fhc1h.lbound(1); j <= fhc1h.ubound(1); ++j) {
	for (int i=fhc1h.lbound(0); i <= fhc1h.ubound(0); ++i) {
//		double fg1 = fgice1(i,j);
//		double fg1g = fgice1_glint2(i,j);
//
//		if ((fg1 > 0) && (fg1 != fg1g)) {
		if (fgice1(i,j) > 0) {
			double val = 1.0d - fgice1_glint2(i,j) / fgice1(i,j);
			if (std::abs(val) < 1e-13) val = 0;
			fhc1h(i,j,1) = val;
		}
	}}

printf("init_landice_com_part2 3\n");
	// ======================= fhc(:,:,hp>1)
	HCIndex hc_index(api->maker->n1());
	std::unique_ptr<giss::VectorSparseMatrix> hp_to_atm(api->maker->hp_to_atm());
	ModelEDomain &domain(*api->domain);

	// Filter this array, and convert to fhc format
	for (auto ii = hp_to_atm->begin(); ii != hp_to_atm->end(); ++ii) {
		int lindex1a[domain.num_local_indices];

		// Input: HP space
		int lindex[domain.num_local_indices];
		int hp1b, i1b;
		hc_index.index_to_ik(ii.col(), i1b, hp1b);
		domain.global_to_local(i1b, lindex);
		if (!domain.in_domain(lindex)) continue;

		// Output: GCM grid
		int i1a = ii.row();
		if (i1a != i1b) {
			fprintf(stderr, "HP2ATM matrix is non-local!\n");
			throw std::exception();
		}

		// Now fill in FHC
		// +1 for C-to-Fortran conversion
		// +1 because lowest HP/HC is reserved for non-model ice
		fhc1h(lindex[0], lindex[1], hp1b+2) +=
			ii.val() * (1.0d - fhc1h(lindex[0], lindex[1],1));
	}
	hp_to_atm.release();

printf("init_landice_com_part2 4\n");
	// ====================== used
	auto used1h(used1h_f.to_blitz());
	used1h = 0;

	for (int j=fhc1h.lbound(1); j <= fhc1h.ubound(1); ++j) {
	for (int i=fhc1h.lbound(0); i <= fhc1h.ubound(0); ++i) {
		// Nothing to do if there's no ice in this grid cell
		if (fgice1(i,j) == 0) continue;

		// Set used for the legacy height point
		// Compute legacy height point for ALL cells with ice
		used1h(i,j,1) = 1;
		// Compute legacy height point just for cells with non-model ice
		// used1h(i,j,1) = (fhc1h(i,j,1) > 0 ? 1 : 0);

		// Min & max height point used for each grid cell
		int mink = std::numeric_limits<int>::max();
		int maxk = std::numeric_limits<int>::min();

		// Loop over HP's (but not the reserved ones) to find
		// range of HP's used on this grid cell.
		for (int k=2; k <= nhp; ++k) {
			if (fhc1h(i,j,k) > 0) {
				mink = std::min(mink, k);
				maxk = std::max(maxk, k);
			}
		}

		// Add a couple of HPs around it!
		mink = std::max(2, mink-2);
		maxk = std::min(nhp, maxk+2);

		// Set everything from mink to maxk as used
		for (int k=mink; k<maxk; ++k) used1h(i,j,k) = 1;
	}}

	// ModelE hack: ModelE disregards used1h, it turns on a height point
	// iff fhc != 0.  So make sure fhc is non-zero everywhere usedhp is set.
	for (int k=fhc1h.lbound(2); k <= fhc1h.ubound(2); ++k) {
	for (int j=fhc1h.lbound(1); j <= fhc1h.ubound(1); ++j) {
	for (int i=fhc1h.lbound(0); i <= fhc1h.ubound(0); ++i) {
		if (used1h(i,j,k) && (fhc1h(i,j,k) == 0)) fhc1h(i,j,k) = 1e-30;
	}}}


printf("END glint2_modele_init_landice_com_part2\n");
}

extern "C"
void glint2_modele_init_hp_to_ices(glint2::modele::glint2_modele *api)
{
printf("BEGIN glint2_modele_init_hp_to_ices\n");
	ModelEDomain &domain(*api->domain);
	HCIndex hc_index(api->maker->n1());

	// ====================== hp_to_ices
	api->hp_to_ices.clear();
	for (auto sheet=api->maker->sheets.begin(); sheet != api->maker->sheets.end(); ++sheet) {

		// Get matrix for HP2ICE
		std::unique_ptr<giss::VectorSparseMatrix> imat(sheet->hp_to_ice());
		if (imat->size() == 0) continue;

		// Convert to GCM coordinates
		std::vector<hp_to_ice_rec> omat;
		omat.reserve(imat->size());
		for (auto ii=imat->begin(); ii != imat->end(); ++ii) {
			// Get index in HP space
			int lindex[domain.num_local_indices];
			int hp1, i1;
			hc_index.index_to_ik(ii.col(), i1, hp1);
			domain.global_to_local(i1, lindex);
			if (!domain.in_domain(lindex)) continue;

			// Write to output matrix
			// +1 for C-to-Fortran conversion
			// +1 because lowest HP/HC is reserved for non-model ice
			omat.push_back(hp_to_ice_rec(
				ii.row(),
				lindex[0], lindex[1], hp1+2,
				ii.val()));
		}

		// Store away
		api->hp_to_ices[sheet->index] = std::move(omat);
	}

printf("END glint2_modele_init_hp_to_ices\n");
}
// -----------------------------------------------------
/** @param hpvals Values on height-points GCM grid for various fields
	the GCM has decided to provide. */
extern "C"
void glint2_modele_couple_to_ice_c(
glint2_modele *api,
int itime,
giss::F90Array<double,3> &smb1h_f,
giss::F90Array<double,3> &seb1h_f,
giss::F90Array<double,3> &tg21h_f)
{
	GCMCoupler_MPI &coupler(*api->gcm_coupler);
int rank = coupler.rank();	// debugging

	std::vector<IceField> fields =
		{IceField::MASS_FLUX, IceField::ENERGY_FLUX, IceField::TG2};
//	std::vector<blitz::Array<double,3>> vals1hp =
//		{smb1h_f.to_blitz(), seb1h_f.to_blitz()};

	auto smb1h(smb1h_f.to_blitz());
	auto seb1h(seb1h_f.to_blitz());
	auto tg21h(tg21h_f.to_blitz());

	// Count total number of elements in the matrices
	// (_l = local to this MPI node)
	int nele_l = 0; //api->maker->ice_matrices_size();
	for (auto ii = api->hp_to_ices.begin(); ii != api->hp_to_ices.end(); ++ii) {
		nele_l += ii->second.size();
	}

	// Allocate buffer for that amount of stuff
	int nfields = fields.size();
	giss::DynArray<SMBMsg> sbuf(SMBMsg::size(nfields), nele_l);

	// Fill it in by doing a sparse multiply...
	// (while translating indices to local coordinates)
	HCIndex hc_index(api->maker->n1());
	int nmsg = 0;
printf("[%d] hp_to_ices.size() = %ld\n", rank, api->hp_to_ices.size());
	for (auto ii = api->hp_to_ices.begin(); ii != api->hp_to_ices.end(); ++ii) {
		int sheetno = ii->first;
		std::vector<hp_to_ice_rec> &mat(ii->second);

printf("[%d] mat[sheetno=%d].size() == %ld\n", rank, sheetno, mat.size());
		// Skip if we have nothing to do for this ice sheet
		if (mat.size() == 0) continue;

		// Do the multiplication
		for (int j=0; j < mat.size(); ++j) {
			hp_to_ice_rec &jj(mat[j]);

			SMBMsg &msg = sbuf[nmsg];
			msg.sheetno = sheetno;
			msg.i2 = jj.row;

			msg[0] = jj.val * smb1h(jj.col_i, jj.col_j, jj.col_k);
			msg[1] = jj.val * seb1h(jj.col_i, jj.col_j, jj.col_k);
			msg[2] = jj.val * tg21h(jj.col_i, jj.col_j, jj.col_k);

//printf("msg = %d (i,j, hc)=(%d %d %d) i2=%d %g %g (%g %g)\n", msg.sheetno, lindex[0], lindex[1], ihc+1, msg.i2, msg[0], msg[1], smb1h(lindex[0], lindex[1], ihc+1), seb1h(lindex[0], lindex[1], ihc+1));

			++nmsg;
		}
	}

	// Sanity check: make sure we haven't overrun our buffer
	if (nmsg != sbuf.size) {
		fprintf(stderr, "Wrong number of items in buffer: %d vs %d expected\n", nmsg, sbuf.size);
		throw std::exception();
	}

	coupler.couple_to_ice(itime, fields, sbuf);
}
