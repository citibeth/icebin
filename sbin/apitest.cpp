#include <glint2/modele_api.hpp>

using namespace glint2;
using namespace glint2::modele;

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

	modele_api *api = modele_api_new(
		maker_fname.c_str(), maker_fname.size(),
		maker_vname.c_str(), maker_vname.size(),

		im, jm,

		// int i0h, int i1h, int j0h, int j1h,
		1, im, 0, jm+1,

		// int i0, int i1, int j0, int j1,
		1, im, 1, jm,

		// int j0s, int j1s,
		2, jm-1,

		// MPI Stuff
		// int comm_f, int root;
		MPI_Comm_c2f(comm), 0);

	int nhc = api->maker->nhc();

	auto fhc1h = blitz::Array<double,3>(
		blitz::Range(1,im),
		blitz::Range(1,jm),
		blitz::Range(1,nhc),
		blitz::fortranArray);
	giss::F90Array<double,3> fhc1h_f(fhc1h);

	auto fgice1 = blitz::Array<double,2>(
		blitz::Range(1,im),
		blitz::Range(1,jm),
		blitz::fortranArray);
	giss::F90Array<double,2> fgice1_f(fgice1);


	modele_api_compute_fhc_c(api, fhc1h_f, fgice1_f);


	modele_api_delete(api);

	// -----------------------------------
		
	printf("Goodbye world from processor %s, rank %d"
		" out of %d processors\n",
		processor_name, world_rank, world_size);

	// Finalize the MPI environment.
	MPI_Finalize();

	return 0;
}
