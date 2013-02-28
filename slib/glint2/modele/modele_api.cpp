#include <glint2/modele/ModelEDomain.hpp>
#include <giss/f90blitz.hpp>

using namespace glint2;
using namespace glint::modele;

struct modele_api {
	std::unique_ptr<MatrixMaker> maker;
	ModelEDomain *domain;
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
extern "C" modele_api_compute_fhc_(modele_api *api,
	giss::F90Array<double, 3> &fhc1h_f,
	giss::F90Array<double, 2> &fgice1_f)
{
	// Reconstruct arrays, using Fortran conventions
	// (smallest stride first, whatever-based indexing it came with)
	auto fhc1h(fhc1h_f.to_blitz());
	auto fgice1(fgice1_f.to_blitz());

	// Zero out fhc1h and fgice1...
	fhc1h = 0;
	fgice1 = 0;
#if 0
	ModelEDomain const &domain(*api->domain);
	for (int j=domain.j0h_f; j <= domain.j1h_f; ++j) {
	for (int i=domain.i0h_f; i <= domain.i1h_f; ++i) {
		fhc1h(i,j) = 0;
		fgice1(i,j) = 0;
	}}
#endif

	// Get the sparse vector values
	std::vector<int> &indices1;	// i1
	std::vector<double> &fhc1h_vals;	// [*nhc]
	std::vector<double> &fgice1_vals;
	compute_fhc2(indices1, fhc1h_vals, fgice1_vals);

	// Put them into the Fortran arrays for this domain
	for (size_t ii=0; ii<indices1.size(); ++ii) {
		int i1 = indices1[ii];

		// Filter out things not in our halo
		if (!domain.in_halo(i1)) continue;

		// Convert to local (ModelE 2-D) indexing convention
		int local[domain.num_local_indices];
		domain.global_to_local(i1, local);

		// Add into Fortran array
		fhc1h(local[0], local[1]) += fhc1h_vals[ii];
		fgice1(local[0], local[1]) += fgice1_vals[ii];
	}
}
// -----------------------------------------------------
// -----------------------------------------------------
// -----------------------------------------------------
