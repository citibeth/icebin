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
	std::string constants_fname = "modele_constants.nc";
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
	int iyear1 = giss::get_att(rparam, "iyear1")->as_int(0);
	printf("dtsrc = %f\n", dtsrc);
	printf("iyear1 = %d\n", iyear1);

	// -----------------------------------------------------
	// Initialize GLINT2 (via a Fortran-ish interface)
	NcVar *time_nc = dnc.get_var("time");
	long cur[4] = {0, 0, 0, 0};
	long counts1 = 1;
	time_nc->set_cur(cur);	// Just use cur[0]
	double time_si;		// itimei
	time_nc->get(&time_si, &counts1);	// Just get one item
	int itimei = (time_si / dtsrc + .5);

	glint2_modele *api = new_glint2_modele();

	// Set up constants from ModelE here...
	NcFile cnc(constants_fname.c_str());
	giss::ConstantSet &gcm_constants(api->gcm_coupler->gcm_constants);
	gcm_constants.read_from_netcdf(cnc, "constants");
	cnc.close();

	glint2_modele_init0(api,
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

		// MPI Stuff
		// int comm_f, int root;
		MPI_Comm_c2f(comm), 0,

		// API Control: write_constants = false
		0);

	glint2_modele_set_start_time(api, iyear1, itimei, dtsrc);

	int nhp = api->maker->nhp(-1);
	printf("desm.cpp: Number Elevation Points (nhp)== %d\n", nhp);

	glint2_modele_init_hp_to_ices(api);

	// ----------------------------------------------------------
	// -------- Test regridding to the ice model

	// Get Time...
//	NcVar *date_nc = dnc.get_var("date");

	// Get Data
	const int nvar = 3;
	NcVar *vars_nc[nvar] = {
		dnc.get_var("impm_lndice"),
		dnc.get_var("imph_lndice"),
		dnc.get_var("tice_lndice")};

	// Get dimensions by querying one variable
	NcVar *var_nc = vars_nc[0];
	long ntime = var_nc->get_dim(0)->size();
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

	blitz::Array<double,3> &impm_c(vars_c[0]);
	blitz::Array<double,3> &imph_c(vars_c[1]);
	blitz::Array<double,3> &tice_c(vars_c[2]);


	// The main loop
	for (long time_i=0; time_i<ntime-1; ++time_i) {
		// End of this time interval
		cur[1] = time_i+1;
		time_nc->set_cur(cur);	// Just use cur[0]
		double end_time_s;
		time_nc->get(&end_time_s, counts);

		// Beginning of this time interval.
		cur[0] = time_i;
		time_nc->set_cur(cur);	// Just use cur[0]
		double begin_time_s;
		time_nc->get(&begin_time_s, counts);

		// Copuling time interval (in seconds)
		double dt_s = end_time_s - begin_time_s;

		// Leave cur[] set at the beginning of the time interval, for reading variable below

		// itime: ModelE calls GLINT2 at the END of a time interval
		int itime = (end_time_s / dtsrc * .5);

		printf("**** itime=%d, end_time_s=%f\n", itime, end_time_s);

		// Read the variables
		for (int vi=0; vi<nvar; ++vi) {
			NcVar *var_nc = vars_nc[vi];
			blitz::Array<double,3> &var_c = vars_c[vi];
			blitz::Array<double,3> &var_f = vars_f[vi];

			// Read the variable (over the time interval)
			var_nc->set_cur(cur);
			var_nc->get(var_c.data(), counts);

			// Zero out nans
			for (int i=0; i<var_c.extent(0); ++i) {
			for (int j=0; j<var_c.extent(1); ++j) {
			for (int k=0; k<var_c.extent(2); ++k) {
				double &var_ijk = var_c(i,j,k);
				if (fabs(var_ijk) >= 1e10) var_ijk = 0;
			}}}
				
			// HACK: Clear things not in our domain on this MPI node
			for (int j=1; j<j0; ++j)
				var_f(blitz::Range::all(), j, blitz::Range::all()) = 0;
			for (int j=j1+1; j <= jm; ++j)
				var_f(blitz::Range::all(), j, blitz::Range::all()) = 0;

		}

		// ModelE writes out NetCDF files in units X/s.
		// But internally, impm_lndice/etc. are stored in just X.
		// So we multiply by dt_s to get back to the internal ModelE units.
		impm_c *= dt_s;
		imph_c *= dt_s;

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
