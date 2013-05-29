#include <mpi.h>		// Must come first for Intel MPI
#include <boost/function.hpp>
#include <netcdfcpp.h>
#include <glint2/modele/glint2_modele.hpp>
#include <giss/f90blitz.hpp>

using namespace glint2;
using namespace glint2::modele;

// --------------------------------------------------------
int main(int argc, char **argv)
{

	// Initialize the MPI environment
	MPI_Init(&argc, &argv);

	// Get the communicator
	MPI_Comm comm = MPI_COMM_WORLD;
 
	// Get the number of processes
	int world_size;
	MPI_Comm_size(comm, &world_size);
 
	// Get the rank of the process
	int world_rank;
	MPI_Comm_rank(comm, &world_rank);
 
	// Get the name of the processor
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
 
	// Print off a hello world message
	printf("Hello world from processor %s, rank %d"
		" out of %d processors\n",
		processor_name, world_rank, world_size);
 
	// -----------------------------------
	std::string maker_fname = "mmx.nc";
	std::string maker_vname = "m";

	int im = 144;
	int jm = 90;

	int j0, j1;
	int j0s, j1s;
	switch(world_rank) {
		case 0 :
			j0 = 1;
			j1 = 79;
			j0s = 2;
			j1s = 79;
		break;
		case 1 :
			j0 = 80;
			j1 = jm;
			j0s = 80;
			j1s = jm-1;
		break;
	}

#if 0
	j0 = 1;
	j1 = jm;
	j0s = 2;
	j1s = jm-1;
#endif




	glint2_modele *api = glint2_modele_new(
		maker_fname.c_str(), maker_fname.size(),
		maker_vname.c_str(), maker_vname.size(),

		im, jm,

		// int i0h, int i1h, int j0h, int j1h,
		1, im, j0-1, j1+1,

		// int i0, int i1, int j0, int j1,
		1, im, j0, j1,

		// int j0s, int j1s,
		// (Start and end of domain exclusive of poles
		j0s, j1s,
//		2, jm-1,

		// MPI Stuff
		// int comm_f, int root;
		MPI_Comm_c2f(comm), 0);

	int nhp = api->maker->nhp();

	auto used1h = blitz::Array<int,3>(
		blitz::Range(1,im),
		blitz::Range(1,jm),
		blitz::Range(1,nhp),
		blitz::fortranArray);
	giss::F90Array<int,3> used1h_f(used1h);

	auto fhc1h = blitz::Array<double,3>(
		blitz::Range(1,im),
		blitz::Range(1,jm),
		blitz::Range(1,nhp),
		blitz::fortranArray);
	giss::F90Array<double,3> fhc1h_f(fhc1h);

#define FRAC_VAR(name) \
	auto name = blitz::Array<double,2>( \
		blitz::Range(1,im), \
		blitz::Range(1,jm), \
		blitz::fortranArray); \
	giss::F90Array<double,2> name##_f(name)

	FRAC_VAR(fgice1_glint2);
	FRAC_VAR(fgice1);
	FRAC_VAR(fgrnd1);
	FRAC_VAR(focean1);
	FRAC_VAR(flake1);

	{
		NcFile nc("Z2HX2fromZ1QX1N_hc.nc");
		long counts[2] = {jm, im};
		nc.get_var("fgice")->get(fgice1.data(), counts);
		nc.get_var("fgrnd")->get(fgrnd1.data(), counts);
		nc.get_var("focean")->get(focean1.data(), counts);
		nc.get_var("flake")->get(flake1.data(), counts);
		nc.close();
	}

	for (int j=1; j<j0; ++j) {
		fgice1_glint2(blitz::Range::all(), j) = 0;
		fgice1(blitz::Range::all(), j) = 0;
		fgrnd1(blitz::Range::all(), j) = 0;
		focean1(blitz::Range::all(), j) = 0;
		flake1(blitz::Range::all(), j) = 0;
	}
	for (int j=j1+1; j <= jm; ++j) {
		fgice1_glint2(blitz::Range::all(), j) = 0;
		fgice1(blitz::Range::all(), j) = 0;
		fgrnd1(blitz::Range::all(), j) = 0;
		focean1(blitz::Range::all(), j) = 0;
		flake1(blitz::Range::all(), j) = 0;
	}


	glint2_modele_compute_fgice_c(api, 1 /*true*/,
		fgice1_glint2_f,
		fgice1_f, fgrnd1_f, focean1_f, flake1_f);

	// ----------------------------------------------------------
	// Try the HP-to-HC matrix
	blitz::Array<double,3> elevhc(
		blitz::Range(1,im),
		blitz::Range(1,jm),
		blitz::Range(1,nhp),
		blitz::fortranArray);
	giss::F90Array<double, 3> elevhc_f(elevhc);
	blitz::Array<double,3> fhp_approx1h(
		blitz::Range(1,im),
		blitz::Range(1,jm),
		blitz::Range(1,nhp),
		blitz::fortranArray);
	giss::F90Array<double, 3> fhp_approx1h_f(fhp_approx1h);

	blitz::Array<double,2> zatmo1(
		blitz::Range(1,im),
		blitz::Range(1,jm),
		blitz::fortranArray);
	giss::F90Array<double, 2> zatmo1_f(zatmo1);

	zatmo1 = 0;

	glint2_modele_init_landice_com_c(api,
		zatmo1_f, 1.0, fgice1_glint2_f, fgice1_f,
		used1h_f, fhc1h_f, elevhc_f);

	// ----------------------------------------------------------
	// Save it to a netCDF file so we can tell if it's correct
	char fname[50];
	sprintf(fname, "fhc-%02d.nc", world_rank);
	NcFile ncout(fname, NcFile::Replace);

	// Define variables
	std::vector<boost::function<void ()>> fns;
	NcDim *im_dim = ncout.add_dim("im", im);
	NcDim *jm_dim = ncout.add_dim("jm", jm);
	NcDim *nhp_dim = ncout.add_dim("nhp", nhp);

	auto used1h_c(used1h.transpose(1,0));	// Re-order dimensions for netCDF standard
	fns.push_back(giss::netcdf_define(ncout, "used1h", used1h_c, {nhp_dim, jm_dim, im_dim}));

	auto fhc1h_c(fhc1h.transpose(2,1,0));	// Re-order dimensions for netCDF standard
	fns.push_back(giss::netcdf_define(ncout, "fhc1h", fhc1h_c, {nhp_dim, jm_dim, im_dim}));

	auto fgice1_c(fgice1.transpose(1,0));
	fns.push_back(giss::netcdf_define(ncout, "fgice1", fgice1_c, {jm_dim, im_dim}));

	auto fgrnd1_c(fgrnd1.transpose(1,0));
	fns.push_back(giss::netcdf_define(ncout, "fgrnd1", fgrnd1_c, {jm_dim, im_dim}));

	auto focean1_c(focean1.transpose(1,0));
	fns.push_back(giss::netcdf_define(ncout, "focean1", focean1_c, {jm_dim, im_dim}));

	auto flake1_c(flake1.transpose(1,0));
	fns.push_back(giss::netcdf_define(ncout, "flake1", flake1_c, {jm_dim, im_dim}));

	// ----------------------------------------------------------
	// -------- Test regridding to the ice model
//	NcFile ijhc_nc("JUL1952.ijhcncar225.nc");
	NcFile ijhc_nc("JUL1950.ijhchc1k225.nc");

	auto impm_c(giss::read_blitz<double,3>(ijhc_nc, "impm_lndice"));
	auto impm(giss::c_to_f(impm_c));	// Fortran-style array
	giss::F90Array<double,3> impm_f(impm);

	auto imph_c(giss::read_blitz<double,3>(ijhc_nc, "imph_lndice"));
	auto imph(giss::c_to_f(imph_c));	// Fortran-style array
	giss::F90Array<double,3> imph_f(imph);


	for (int i=0; i<impm_c.extent(0); ++i) {
	for (int j=0; j<impm_c.extent(1); ++j) {
	for (int k=0; k<impm_c.extent(2); ++k) {
		double &impm_ijk = impm_c(i,j,k);
		if (fabs(impm_ijk) >= 1e10) impm_ijk = 0;
		double &imph_ijk = imph_c(i,j,k);
		if (fabs(imph_ijk) >= 1e10) imph_ijk = 0;
	}}}



	ijhc_nc.close();

	// HACK: Clear things not in our domain
	for (int j=1; j<j0; ++j) {
		impm(blitz::Range::all(), j, blitz::Range::all()) = 0;
		imph(blitz::Range::all(), j, blitz::Range::all()) = 0;
	}
	for (int j=j1+1; j <= jm; ++j) {
		impm(blitz::Range::all(), j, blitz::Range::all()) = 0;
		imph(blitz::Range::all(), j, blitz::Range::all()) = 0;
	}

//	glint2_modele_couple_to_ice(api, impm_f, imph_f);
	glint2_modele_couple_to_ice_c(api, impm_f);
	
	ijhc_nc.close();



	// --------------------------------------------------------

	// --------------------------------------------------------
	// Write it out
	for (auto ii = fns.begin(); ii != fns.end(); ++ii) (*ii)();
printf("AA\n");

	// Finish up
	glint2_modele_delete(api);
	ncout.close();

	// -----------------------------------
		
	printf("Goodbye world from processor %s, rank %d"
		" out of %d processors\n",
		processor_name, world_rank, world_size);

	// Finalize the MPI environment.
	MPI_Finalize();

	return 0;
}
