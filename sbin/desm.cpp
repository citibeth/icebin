/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <mpi.h>		// Must come first for Intel MPI
#include <boost/function.hpp>
#include <netcdfcpp.h>
#include <glint2/modele/glint2_modele.hpp>
#include <giss/f90blitz.hpp>
#include <boost/filesystem.hpp>

using namespace glint2;
using namespace glint2::modele;

// !@param shi heat capacity of pure ice (at 0 C) (2060 J/kg C)
const double SHI  = 2060.;
//@param lhm   latent heat of melt at 0 C (334590 J/kg)
const double LHM = 3.34e5;


// --------------------------------------------------------
int main(int argc, char **argv)
{
	std::string glint2_config_fname = "modele_ll_g2x2_5-searise_g20-40-PISM.nc";
//	std::string glint2_config_fname = "modele_ll_g2x2_5-searise_g20-40-DISMAL.nc";
	std::string glint2_config_vname = "m";
	std::string desm_fname = "desm_in.nc";

	// Initialize the MPI environment
	MPI_Init(&argc, &argv);

	// Get the communicator
	MPI_Comm comm = MPI_COMM_WORLD;
 
	// Get the number of processes
	int world_size;
	MPI_Comm_size(comm, &world_size);
 	if (world_size > 2) {
		printf("This test program not designed to run with -np >2\n");
		return -1;
	}

	// Get the rank of the process
	int world_rank;
	MPI_Comm_rank(comm, &world_rank);
 
	// Get the name of the processor
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
 
	// Print a hello world message
	printf("Hello world from processor %s, rank %d"
		" out of %d processors\n",
		processor_name, world_rank, world_size);
 
	// -----------------------------------
	// Set up domain decomposition

	int im = 144;
	int jm = 90;

	int j0, j1;
	int j0s, j1s;

	if (world_size == 1) {
		j0 = 1;
		j1 = jm;
		j0s = 2;
		j1s = jm-1;
	} else {
		switch(world_rank) {
			// Domains are carefully chosen so Greenland is split between them.
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
	}

	// -------------------------------------------------
	// Open desm input file and read dtsrc (and other parameters)
	printf("desm_fname = %s\n", desm_fname.c_str());
	NcFile dnc(desm_fname.c_str());

	// Get parameters
	NcVar *rparam = dnc.get_var("rparam");
	double dtsrc = giss::get_att(rparam, "dtsrc")->as_double(0);
	printf("dtsrc = %f\n", dtsrc);

	// -----------------------------------------------------
	// Initialize GLINT2 (via a Fortran-ish interface)
	glint2_modele *api = glint2_modele_new(
		glint2_config_fname.c_str(), glint2_config_fname.size(),
		glint2_config_vname.c_str(), glint2_config_vname.size(),

		im, jm,

		// int i0h, int i1h, int j0h, int j1h,
		1, im, j0-1, j1+1,

		// int i0, int i1, int j0, int j1,
		1, im, j0, j1,

		// int j0s, int j1s,
		// (Start and end of domain exclusive of poles
		j0s, j1s,

		// iyear1
		1950,

		// dtsrc  (see MODEL_COM.f)
		dtsrc,

		// MPI Stuff
		// int comm_f, int root;
		MPI_Comm_c2f(comm), 0,

		// Constants
		// (which come directly from ModelE, because this is called from ModelE)
		LHM, SHI);

	int nhp = api->maker->nhp(-1);
	printf("api->maker->nhp(-1) == %d\n", nhp);

	glint2_modele_init_hp_to_ices(api);

	// ----------------------------------------------------------
	// -------- Test regridding to the ice model

	// Get Time...
	NcVar *itime_nc = dnc.get_var("itime");
	NcVar *date_nc = dnc.get_var("date");

	// Get Data
	const int nvar = 3;
	NcVar *vars_nc[nvar] = {
		dnc.get_var("impm_lndice"),
		dnc.get_var("imph_lndice"),
		dnc.get_var("tice_lndice")};

	// Get dimensions by querying one variable
	NcVar *var_nc = vars_nc[0];
	long ntime = var_nc->get_dim(0)->size();
	long cur[4] = {0, 0, 0, 0};
	long counts[4] = {1,				// time
		var_nc->get_dim(1)->size(),		// nhc
		var_nc->get_dim(2)->size(),		// jm
		var_nc->get_dim(3)->size()};		// im

	// Allocate arrays (buffers) for one timestep
	blitz::Array<double, 3> vars_c[nvar];		// C ordering of dimensions, 0-based
	blitz::Array<double, 3> vars_f[nvar];		// Fortran ordering of dimensions, 1-based
	giss::F90Array<double,3> vars_ff[nvar];		// Dope vector
	for (int i=0; i<nvar; ++i) {
		vars_c[i].reference(
			blitz::Array<double, 3>(counts[1], counts[2], counts[3]));
		vars_f[i].reference(giss::c_to_f(vars_c[i]));
		vars_ff[i] = giss::F90Array<double,3>(vars_f[i]);
	}

	giss::F90Array<double,3> &impm_ff(vars_ff[0]);
	giss::F90Array<double,3> &imph_ff(vars_ff[1]);
	giss::F90Array<double,3> &tice_ff(vars_ff[2]);

	// The main loop
	for (long time=0; time<ntime; ++time) {
		// Read itime
		int itime;
		itime_nc->set_cur(cur);	// Just use cur[0]
		itime_nc->get(&itime, counts);	// Just get one item

		double time_s = itime * dtsrc;
		printf("**** itime=%d, time_s=%f\n", itime, time_s);

		// Read the variables
		cur[0] = time;
		for (int vi=0; vi<nvar; ++vi) {
			NcVar *var_nc = vars_nc[vi];
			blitz::Array<double,3> &var_c = vars_c[vi];
			blitz::Array<double,3> &var_f = vars_f[vi];

			// Read the variable
			var_nc->set_cur(cur);
			var_nc->get(var_c.data(), counts);

			// Zero out nans
			for (int i=0; i<var_c.extent(0); ++i) {
			for (int j=0; j<var_c.extent(1); ++j) {
			for (int k=0; k<var_c.extent(2); ++k) {
				double &var_ijk = var_c(i,j,k);
				if (fabs(var_ijk) >= 1e10) var_ijk = 0;
			}}}
				
			// HACK: Clear things not in our domain
			for (int j=1; j<j0; ++j)
				var_f(blitz::Range::all(), j, blitz::Range::all()) = 0;
			for (int j=j1+1; j <= jm; ++j)
				var_f(blitz::Range::all(), j, blitz::Range::all()) = 0;

		}

		// Run the coupling step
		glint2_modele_couple_to_ice_c(api, itime, impm_ff, imph_ff, tice_ff);
	}

	// Finish up
	glint2_modele_delete(api);

	// -----------------------------------
		
	printf("Goodbye world from processor %s, rank %d"
		" out of %d processors\n",
		processor_name, world_rank, world_size);

	// Finalize the MPI environment.
	MPI_Finalize();

	return 0;
}
