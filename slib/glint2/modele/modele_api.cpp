#include <glint2/modele/ModelEDomain.hpp>
#include <giss/f90blitz.hpp>

using namespace glint2;
using namespace glint::modele;

struct modele_api {
	std::unique_ptr<MatrixMaker> maker;
	ModelEDomain *domain;	// Points to domain owned by maker
	std::unique_ptr<VectorSparseMatrix> hp_to_hc;
}

/** @param spec 
extern "C" modele_api_new_(
	char *maker_fname_f, const int &maker_fname_len,
	char *maker_vname_f, const int &maker_vname_len,

	// Info about the global grid
	const int &im, const int &jm,

	// Info about the local grid (C-style indices)
	const int &i0h, const int &i1h, const int &j0h, const int &j1h,
	const int &i0, const int &i1, const int &j0, const int &j1,
	const int &j0s, const int &j1s)
{
	// Convert Fortran arguments
	std::string maker_fname(maker_fname_f, maker_fname_len);
	std::string maker_vname(maker_vname_f, maker_vname_len);

	// Allocate our return variable
	std::unique_ptr<modele_api> api(new model_api());

	// Set up the domain
	std::unique_ptr<GridDomain> mdomain(
	new ModelEDomain(
		i0h, i1h, j0h, j1h,
		i0, i1, j0, j1,
		j0s, j1s));
	api->domain = (ModelEDomain *)mdomain.get();

	// Load the MatrixMaker	(filtering by our domain, of course)
	api->maker.reset(new MatrixMaker(std::move(mdomain)));
	NcFile nc(maker_fname, NcFile::ReadOnly);
	api->maker->read_from_netcdf(nc, maker_vname);
	nc.close();

	// No exception, we can release our pointer back to Fortran
	return api.release();
}
// -----------------------------------------------------
extern "C" modele_api_delete_(modele_api *api)
{
	delete api;
}
// -----------------------------------------------------
extern "C"
void modele_api_compute_fhc_c_(modele_api *api,
	giss::F90Array<double, 3> &fhc1h_f,
	giss::F90Array<double, 2> &fgice1_f)
{
	ModelEDomain &domain(*api->domain);

	// Reconstruct arrays, using Fortran conventions
	// (smallest stride first, whatever-based indexing it came with)
	auto fhc1h(fhc1h_f.to_blitz());
	auto fgice1(fgice1_f.to_blitz());

	// Zero out fhc1h and fgice1...
	fhc1h = 0;
	fgice1 = 0;

	// Get the sparse vector values
	giss::CooVector<std::pair<int,int>,double> fhc1h_s
	giss::CooVector<int,double> fgice1_s;
	api->maker->compute_fhc(fhc1h_s, fgice1_s);

	// Translate the sparse vectors to the ModelE data structures
	for (auto ii = fgice1_s.begin(); ii != fgice1_s.end(); ++ii) {
		int i1 = ii->first;

		// Filter out things not in our domain
		// (we'll get the answer for our halo via a halo update)
		if (!domain.in_domain(i1)) continue;

		// Convert to local (ModelE 2-D) indexing convention
		int local[domain.num_local_indices];
		domain.global_to_local(i1, local);

		// Store it away
		// (we've eliminated duplicates, so += isn't needed, but doesn't hurt either)
		fgice1(local[0], local[1]) += ii->second;
	}

	for (auto ii = fhc1h_s.begin(); ii != fhc1h_s.end(); ++ii) {
		int i1 = ii->first.first;
		int hc = ii->first.second;		// zero-based
		double val = ii->second;

		// Filter out things not in our domain
		// (we'll get the answer for our halo via a halo update)
		if (!domain.in_domain(i1)) continue;

		// Convert to local (ModelE 2-D + height class) indexing convention
		int local[domain.num_local_indices];
		domain.global_to_local(i1, local);

		// Store it away
		// (we've eliminated duplicates, so += isn't needed, but doesn't hurt either)
		int hc_f = hc + 1;		// convert zero-based to 1-based arrays
		fhc1h(local[0], local[1], hc_f) += val;
	}
}
// -----------------------------------------------------
/** Call this to figure out how to dimension arrays.
@return Number of elements in the sparse matrix */
extern "C"
int modele_api_hp_to_hc_part1_(modele_api *api)
{
	auto mat(api->maker->hp_to_hc());
	api->hp_to_hc = filter_matrix(*api->domain, *api->domain, *mat);
	return api->hp_to_hc->size();
}
// -----------------------------------------------------
/** Call this after rows, cols and vals have been dimensioned. */
extern "C"
void modele_api_hp_to_hc_part2_(modele_api *api,
	giss::F90Array<double, 3> &rows_i_f,
	giss::F90Array<double, 3> &rows_j_f,
	giss::F90Array<double, 3> &cols_i_f,
	giss::F90Array<double, 3> &cols_j_f,
	giss::F90Array<double, 3> &vals_f)
{

	// Array bounds checking not needed, it's done
	// in lower-level subroutines that we call.

	// Copy the rows while translating
	std::vector<blitz::Array<int,1>> lrows = {rows_i_f.to_blitz(), rows_j_f.to_blitz()};
	global_to_local(*api->domain, *api->hp_to_hc.rows(), lrows);

	// Copy the cols while translating
	std::vector<blitz::Array<int,1>> lcols = {cols_i_f.to_blitz(), cols_j_f.to_blitz()};
	global_to_local(*api->domain, *api->hp_to_hc.cols(), lcols);

	// Copy the values, just a simple vector copy
	auto vals(vals_f.to_blitz());
	vals = vector_to_blitz(api->hp_to_hc.vals());

	// Free temporary storage
	api->hp_to_hc.reset();
}
// -----------------------------------------------------
