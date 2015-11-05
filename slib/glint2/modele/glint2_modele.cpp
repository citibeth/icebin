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

#include <mpi.h>	// For Intel MPI, mpi.h must be included before stdio.h
#include <netcdfcpp.h>
#include <giss/mpi.hpp>
#include <giss/blitz.hpp>
#include <giss/f90blitz.hpp>
#include <glint2/HCIndex.hpp>
#include <glint2/modele/glint2_modele.hpp>
#include <glint2/modele/GCMCoupler_ModelE.hpp>
//#include <glint2/IceModel_TConv.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <glint2/contracts/contracts.hpp>
#include <giss/exit.hpp>

using namespace glint2;
using namespace glint2::modele;


// ================================================================

struct ModelEMsg {
	int i, j, k;	// Indices into ModelE
	double vals[1];		// Always at least one val; but this could be extended, based on # of inputs

	double &operator[](int i) { return *(vals + i); }

	/** @return size of the struct, given a certain number of values */
	static size_t size(int nfields)
		{ return sizeof(ModelEMsg) + (nfields-1) * sizeof(double); }

	static MPI_Datatype new_MPI_struct(int nfields);

	/** for use with qsort */
//	static int compar(void const * a, void const * b);

};

MPI_Datatype ModelEMsg::new_MPI_struct(int nfields)
{
	int nele = 3 + nfields;
	int blocklengths[] = {1, 1, 1, nfields};
	MPI_Aint displacements[] = {offsetof(ModelEMsg,i), offsetof(ModelEMsg,j), offsetof(ModelEMsg,k), offsetof(ModelEMsg, vals)};
	MPI_Datatype types[] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE};
	MPI_Datatype ret;
	MPI_Type_create_struct(4, blocklengths, displacements, types, &ret);
	MPI_Type_commit(&ret);
	return ret;
}
// -----------------------------------------------------

/** LOCAL FUNCTION: Set up a netCDF file, ready to store timeseries
variables in by Glint2.  This is used for both gcm_input and gcm_ouput
files. */
extern "C"
void init_ncfile(glint2_modele *api,
giss::CouplingContract const &contract,
std::string const &fname)
{
	printf("BEGIN init_ncfile(%s)\n", fname.c_str());
	GCMParams const &gcm_params(api->gcm_coupler.gcm_params);

	// Set up NetCDF file to store GCM output as we received them (modele_out.nc)
	int nhp_gcm = glint2_modele_nhp_gcm(api);
	NcFile ncout(fname.c_str(), NcFile::Replace);
	NcDim *im_dim = ncout.add_dim("im", api->domain->im);
	NcDim *jm_dim = ncout.add_dim("jm", api->domain->jm);
	NcDim *nhp_dim = ncout.add_dim("nhp", nhp_gcm);
	NcDim *one_dim = ncout.add_dim("one", 1);
	NcDim *time_dim = ncout.add_dim("time");		// No dimsize --> unlimited

	const NcDim *dims_elevation[4]{time_dim, nhp_dim, jm_dim, im_dim};
	const NcDim *dims_atmosphere[3]{time_dim, jm_dim, im_dim};

	const NcDim *dims_b[1]{one_dim};
	NcVar *grid_var = ncout.add_var("grid", giss::get_nc_type<double>(), 1, dims_b);
	grid_var->add_att("file", api->gcm_coupler.fname.c_str());
	grid_var->add_att("variable", api->gcm_coupler.vname.c_str());

	NcVar *time0_var = ncout.add_var("time0", giss::get_nc_type<double>(), 1, dims_b);
	time0_var->add_att("units", gcm_params.time_units.c_str());
	time0_var->add_att("calendar", "365_day");
	time0_var->add_att("axis", "T");
	time0_var->add_att("long_name", "Simulation start time");

	NcVar *time_var = ncout.add_var("time", giss::get_nc_type<double>(), 1, dims_atmosphere);
	time_var->add_att("units", gcm_params.time_units.c_str());
	time_var->add_att("calendar", "365_day");
	time_var->add_att("axis", "T");
	time_var->add_att("long_name", "Coupling times");

	//giss::CouplingContract &gcm_outputs(api->gcm_coupler.gcm_outputs);

	for (unsigned int i=0; i < contract.size_nounit(); ++i) {

		NcVar *nc_var = 0;
		giss::CoupledField const &cf(contract.field(i));
printf("Creating NC variable for %s (%s)\n", cf.name.c_str(), contracts::to_str(cf.flags).c_str());
		switch(cf.flags & contracts::GRID_BITS) {
			case contracts::ATMOSPHERE :
				nc_var = ncout.add_var(cf.name.c_str(),
					giss::get_nc_type<double>(), 3, dims_atmosphere);
			break;
			case contracts::ELEVATION :
				nc_var = ncout.add_var(cf.name.c_str(),
				giss::get_nc_type<double>(), 4, dims_elevation);
			break;
			default:
				fprintf(stderr, "init_ncfile() unable to handle grid type %s for field %s\n", contracts::to_str(cf.flags).c_str(), cf.name.c_str());
				giss::exit(1);
		}

		auto comment(boost::format(
			"%s[t,...] holds the mean from time[t-1] to time[t].  See time0[0] if t=0.")
			% contract.name(i));
		nc_var->add_att("comment", comment.str().c_str());

		std::string const &description(contract.field(i).get_description());
		if (description != "") nc_var->add_att("long_name", description.c_str());

		std::string const &units(contract.field(i).get_units());
		if (units != "") nc_var->add_att("units", units.c_str());
	}

	// Put initial time in it...
	long cur[1]{0};
	long counts[1]{1};
	time0_var->set_cur(cur);
	time0_var->put(&gcm_params.time_start_s, counts);

	ncout.close();
	printf("END init_ncfile(%s)\n", fname.c_str());

}
// -----------------------------------------------------

/** LOCAL FUNCTION: Save the ice model output, after it's been
transformed by Glint2 (now it's GCM input). */
static void save_gcm_inputs(
glint2_modele *api,
double time_s,
giss::F90Array<double,3> &gcm_inputs_d_f)	// Fortran array
{
	printf("BEGIN save_gcm_inputs(%s)\n", api->gcm_coupler.gcm_in_file.c_str());

	// Fortran-style array: i,j,ihp, indexing starts at 1
	blitz::Array<double,3> gcm_inputs_d(gcm_inputs_d_f.to_blitz());

	// Get dimensions of full domain
	int nhp_gcm = glint2_modele_nhp_gcm(api);
	GCMCoupler &coupler(api->gcm_coupler);
	giss::CouplingContract const &contract(coupler.gcm_inputs);

	// Open output netCDF file
	NcFile ncout(api->gcm_coupler.gcm_in_file.c_str(), NcFile::Write);	// Read/Write
	NcDim *time_dim = ncout.get_dim("time");
	NcDim *nhp_dim = ncout.get_dim("nhp");
	NcDim *jm_dim = ncout.get_dim("jm");
	NcDim *im_dim = ncout.get_dim("im");

	long cur_ijhc[4]{time_dim->size(),0,0,0};		// time, nhp_gcm, jm, im
	long counts_ijhc[4]{1, nhp_dim->size(), jm_dim->size(), im_dim->size()};

	long cur_ij[4]{time_dim->size(),0,0};		// time, nhp, jm, im
	long counts_ij[4]{1, jm_dim->size(), im_dim->size()};

	NcVar *time_var = ncout.get_var("time");
	time_var->set_cur(cur_ijhc);	// Could be cur_ij, either way is fine
	time_var->put(&time_s, counts_ijhc);

	// Write arrays to it
	int base_index = 0;
	for (unsigned int i=0; i < contract.size_nounit(); ++i) {
		double const *array_base = &gcm_inputs_d(1,1,base_index);	// i,j,ihp

		NcVar *nc_var = ncout.get_var(contract.name(i).c_str());

		switch(contract.field(i).flags & contracts::GRID_BITS) {
			case contracts::ATMOSPHERE :
				nc_var->set_cur(cur_ij);
				nc_var->put(array_base, counts_ij);
				base_index += 1;
			break;
			case contracts::ELEVATION :
				nc_var->set_cur(cur_ijhc);
				nc_var->put(array_base, counts_ijhc);
				base_index += nhp_gcm;
			break;
			default: ;
		}
	}

	ncout.close();
	printf("END save_gcm_inputs(%s)\n", api->gcm_coupler.gcm_in_file.c_str());
}
// -----------------------------------------------------

/** LOCAL function: Log the GCM output, which is to be passed to the
ice model.
@param hpvals Values on height-points GCM grid for various fields
	the GCM has decided to provide. */
static void save_gcm_outputs(
glint2_modele *api,
double time_s,
std::vector<std::unique_ptr<blitz::Array<double,3>>> &inputs)
{

	// Get dimensions of full domain
	int nhp_gcm = glint2_modele_nhp_gcm(api);
	ModelEDomain const *domain(&*api->domain);

	int const rank = api->gcm_coupler.rank();	// MPI rank; debugging

	printf("[%d] BEGIN save_gcm_outputs(time_s=%f)\n", rank, time_s);

	GCMCoupler &coupler(api->gcm_coupler);

	giss::CouplingContract &gcm_outputs(coupler.gcm_outputs);

	// Count total number of elements in the inputs (for this MPI domain)
	auto &input0(inputs[0]);
	int nele_l =
		(input0->ubound(2) - input0->lbound(2) + 1) *
		(domain->j1_f - domain->j0_f + 1) *
		(domain->i1_f - domain->i0_f + 1);

	// Find the max. number of fields (for input) used for any ice sheet.
	// This will determine the size of our MPI messages.
	int nfields = gcm_outputs.size_nounit();

	// Allocate buffer for that amount of stuff
	giss::DynArray<ModelEMsg> sbuf(ModelEMsg::size(nfields), nele_l);

	// Fill it in....
	int nmsg = 0;
	for (int k=input0->lbound(2); k<=input0->ubound(2); ++k)		// nhp_gcm
	for (int j=domain->j0_f; j <= domain->j1_f; ++j)
	for (int i=domain->i0_f; i <= domain->i1_f; ++i) {
		ModelEMsg &msg = sbuf[nmsg];
		msg.i = i;
		msg.j = j;
		msg.k = k;
		for (unsigned int l=0; l<nfields; ++l) msg[l] = (*inputs[l])(i,j,k);
		++nmsg;
	}

	// Sanity check: make sure we haven't overrun our buffer
	if (nmsg != sbuf.size) {
		fprintf(stderr, "Wrong number of items in buffer: %d vs %d expected\n", nmsg, sbuf.size);
		giss::exit(1);
	}

	// Gather it to root
	GCMParams const &gcm_params(api->gcm_coupler.gcm_params);
	std::unique_ptr<giss::DynArray<ModelEMsg>> rbuf = giss::gather_msg_array(
		gcm_params.gcm_comm, gcm_params.gcm_root, sbuf, nfields, nmsg, 0);

	// Process the gathered data
	if (rank == gcm_params.gcm_root) {

		// Allocate ijk arrays
		std::vector<blitz::Array<double,3>> outputs;
		for (unsigned int i=0; i<nfields; ++i) {
			outputs.push_back(blitz::Array<double,3>(nhp_gcm, domain->jm, domain->im));
		}

		// Turn messages into ijk arrays
		for (auto msg=rbuf->begin(); msg != rbuf->end(); ++msg) {
			for (unsigned int i=0; i<nfields; ++i) {
				int mi = msg->i - 1;
				int mj = msg->j - 1;
				int mk = msg->k - 1;		// Convert Fortran --> C indexing
				outputs[i](mk, mj, mi) = (*msg)[i];
			}
		}

		// Write the arrays to a file
		NcFile ncout(api->gcm_coupler.gcm_out_file.c_str(), NcFile::Write);	// Read/Write
		NcDim *time_dim = ncout.get_dim("time");
		NcDim *nhp_dim = ncout.get_dim("nhp");
		NcDim *jm_dim = ncout.get_dim("jm");
		NcDim *im_dim = ncout.get_dim("im");

		long cur[4]{time_dim->size(),0,0,0};		// time, nhp_gcm, jm, im
		long counts[4]{1, nhp_dim->size(), jm_dim->size(), im_dim->size()};

		NcVar *time_var = ncout.get_var("time");
		time_var->set_cur(cur);
		time_var->put(&time_s, counts);

		for (int i=0; i<nfields; ++i) {
			NcVar *nc_var = ncout.get_var(gcm_outputs.name(i).c_str());
			nc_var->set_cur(cur);
			nc_var->put(outputs[i].data(), counts);
		}

		ncout.close();
	}
	printf("[%d] END save_gcm_outputs(time_s=%f)\n", rank, time_s);
}

// ================================================================
extern "C" glint2_modele *new_glint2_modele_c()
{
	std::unique_ptr<glint2_modele> api(new glint2_modele());

	// No exception was thrown... we can release our pointer back to Fortran
	glint2_modele *ret = api.release();

//	int const rank = ret->gcm_coupler.rank();	// MPI rank; debugging
//	printf("[%d] Allocated glint2_modele api struct: %p\n", rank, ret);
	return ret;
}
// ---------------------------------------------------
/** Set a single constant value in Glint2.  This is a callback, to be called
from ModelE's (Fortran code) constant_set::set_all_constants() */
extern "C" void glint2_modele_set_const(
	glint2::modele::glint2_modele *api,
	char const *name_f, int name_len,
	double val,
	char const *units_f, int units_len,
	char const *description_f, int description_len)
{
//	int const rank = api->gcm_coupler.rank();	// MPI rank; debugging
	giss::ConstantSet &gcm_constants(api->gcm_coupler.gcm_constants);
	gcm_constants.set(
		std::string(name_f, name_len),
		val,
		std::string(units_f, units_len),
		std::string(description_f, description_len));
}

// ---------------------------------------------------
/** Called immediately after glint2_model_set_const...
@param glint2_config_fname_f Name of GLINT2 configuration file */
extern "C" void glint2_modele_init0(
	glint2::modele::glint2_modele *api,
	char const *run_dir_f, int run_dir_len,
	char const *glint2_config_fname_f, int glint2_config_fname_len,
	char const *maker_vname_f, int maker_vname_len,

	// Info about the global grid
	int im, int jm,

	// Info about the local grid (C-style indices)
	int i0h, int i1h, int j0h, int j1h,
	int i0, int i1, int j0, int j1,
	int j0s, int j1s,

	// MPI Stuff
	MPI_Fint comm_f, int root,

	// API  control
	int write_constants)
{
//iyear1=1950;		// Hard-code iyear1 because it hasn't been initialized yet in ModelE

	// Convert Fortran arguments
	std::string run_dir(run_dir_f, run_dir_len);		// The file specified in ModelE -i on the command line (i.e. the processed rundeck)
	std::string glint2_config_fname(glint2_config_fname_f, glint2_config_fname_len);

	printf("BEGIN glint2_modele_init0(%s)\n", glint2_config_fname.c_str());

	std::string maker_vname(maker_vname_f, maker_vname_len);
//glint2_config_fname = "./GLINT2";

	// Parse directory out of glint2_config_fname
	boost::filesystem::path mypath(glint2_config_fname);
	boost::filesystem::path glint2_config_rfname(boost::filesystem::canonical(mypath));
	boost::filesystem::path glint2_config_dir(glint2_config_rfname.parent_path());

	// Set up parmaeters from the GCM to the ice model
	api->gcm_coupler.gcm_params = GCMParams(
		MPI_Comm_f2c(comm_f),
		root,
		glint2_config_dir, run_dir);

	if (write_constants) {
		// Store constants int NetCDF file, so we can desm without ModelE.
		NcFile nc((api->gcm_coupler.gcm_params.run_dir / "modele_constants.nc").c_str(), NcFile::Replace);
		api->gcm_coupler.gcm_constants.netcdf_define(nc, "constants");
		// No NetCDF write phase is required here, constants are all in the define phase.
		nc.close();
	}

#if 1
	// Set up the domain
	std::unique_ptr<GridDomain> mdomain(
		new ModelEDomain(im, jm,
			i0h, i1h, j0h, j1h,
			i0, i1, j0, j1,
			j0s, j1s));
	api->domain = (ModelEDomain *)mdomain.get();

	// ModelE makes symlinks to our real files, which we don't want.

	printf("Opening GLINT2 config file: %s\n", glint2_config_rfname.c_str());

	// Read the coupler, along with ice model proxies
	api->gcm_coupler.read_from_netcdf(glint2_config_rfname.string(), maker_vname, std::move(mdomain));

	// Check bounds on the IceSheets, set up any state, etc.
	// This is done AFTER setup of gcm_coupler because gcm_coupler.read_from_netcdf()
	// might change the IceSheet, in certain cases.
	// (for example, if PISM is used, elev2 and mask2 will be read from related
	// PISM input file, and the version in the GLINT2 file will be ignored)
	api->gcm_coupler.realize();

	// TODO: Test that im and jm are consistent with the grid read.
#endif
	printf("***** END glint2_modele_init0()\n");
}
// -----------------------------------------------------
extern "C" void glint2_modele_delete(glint2_modele *&api)
{
	if (api) delete api;
	api = 0;
}
// -----------------------------------------------------
extern "C"
int glint2_modele_nhp_gcm(glint2_modele const *api)
{
	int nhp_im = api->gcm_coupler.maker->nhp(-1);	// Assume all grid cells have same # EP
	// HP/HC = 1 (Fortran) reserved for legacy "non-model" ice
    // (not part of GLINT2)

	// "nhp_gcm = nhp_im+1" is embedded in the code in this compilation
	// unit.  The first elevation point is skipped in arrays passed from the
	// GCM because that is a non-ice model elevation point.
	return nhp_im + 1;
}
// -----------------------------------------------------
/** @return Number of "copies" of the atmosphere grid must be used
to store the inputs. */
extern "C"
int glint2_modele_gcm_inputs_nhp(glint2_modele *api)
{
	int ihp = api->gcm_inputs_ihp[api->gcm_inputs_ihp.size()-1];
	return ihp;
}
// -----------------------------------------------------
/** @para var_nhp Number of elevation points for this variable.
 (equal to 1 for atmosphere variables, or nhp for elevation-grid variables)
@param return: Start of this variable in the gcm_inputs_local array (Fortran 1-based index) */
extern "C"
int glint2_modele_add_gcm_input(
glint2_modele *api,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
char const *grid_f, int grid_len,
int initial,	// bool
char const *long_name_f, int long_name_len)
{
	std::string field_name(field_name_f, field_name_len);
	std::string units(units_f, units_len);
	std::string grid(grid_f, grid_len);
	std::string long_name(long_name_f, long_name_len);

	int ihp = api->gcm_inputs_ihp[api->gcm_inputs_ihp.size()-1];
	int var_nhp;
	unsigned int flags = 0;
	if (grid == "ATMOSPHERE") {
		var_nhp = 1;
		flags = contracts::ATMOSPHERE;
	} else if (grid == "ELEVATION") {
		var_nhp = glint2_modele_nhp_gcm(api);
		flags = contracts::ELEVATION;
	} else {
		fprintf(stderr, "Unrecognized grid: %s\n", grid.c_str());
		giss::exit(1);
	}

	if (initial) flags |= contracts::INITIAL;

	api->gcm_coupler.gcm_inputs.add_field(field_name, units, flags, long_name);

	api->gcm_inputs_ihp.push_back(ihp + var_nhp);
	int ret = ihp+1;		// Convert to Fortran (1-based) indexing

	printf("glint2_modele_add_gcm_input(%s, %s, %s) --> %d\n", field_name.c_str(), units.c_str(), grid.c_str(), ret);

	return ret;
}
// -----------------------------------------------------
extern "C"
void glint2_modele_set_start_time(glint2_modele *api, int iyear1, int itimei, double dtsrc)
{
	std::cout << "========= GCM Inputs (second time: must be set by now)" << std::endl;
	std::cout << api->gcm_coupler.gcm_inputs << std::endl;

	GCMParams &gcm_params(api->gcm_coupler.gcm_params);

	api->dtsrc = dtsrc;
	double time0_s = itimei * api->dtsrc;

	api->gcm_coupler.set_start_time(
		giss::time::tm(iyear1, 1, 1),
		time0_s);

	api->itime_last = itimei;

	// Finish initialization...
	// -------------------------------------------------
	// Open file to record GCM outputs to Glint2
	if (api->gcm_coupler.gcm_out_file.length() > 0) {
		init_ncfile(api,
			api->gcm_coupler.gcm_outputs,
			api->gcm_coupler.gcm_out_file);
	}

	// Open file to record GCM inputs from Glint2
	if (api->gcm_coupler.gcm_in_file.length() > 0) {
		init_ncfile(api,
			api->gcm_coupler.gcm_inputs,
			api->gcm_coupler.gcm_in_file);
	}
}
// -----------------------------------------------------
extern "C"
void glint2_modele_get_flice_im_c(glint2_modele *api,
	giss::F90Array<double, 2> &flice1_im_f)		// OUT
{

printf("BEGIN glint2_modele_get_flice_im_c()\n");
std::cout << "flice1_im_f: " << flice1_im_f << std::endl;

	ModelEDomain &domain(*api->domain);

	// Reconstruct arrays, using Fortran conventions
	// (smallest stride first, whatever-based indexing it came with)

	// Get the sparse vector values
	giss::VectorSparseVector<int,double> flice1_s;
	api->gcm_coupler.maker->fgice(flice1_s);

	// Translate the sparse vectors to the ModelE data structures
	std::vector<std::tuple<int, int, double>> flice1_vals;
	for (auto ii = flice1_s.begin(); ii != flice1_s.end(); ++ii) {
		int i1 = ii->first;

		// Filter out things not in our domain
		// (we'll get the answer for our halo via a halo update)
		// Convert to local (ModelE 2-D) indexing convention
		int lindex[domain.num_local_indices];
		domain.global_to_local(i1, lindex);
		if (!domain.in_domain(lindex)) continue;

		// Store it away
		// (we've eliminated duplicates, so += isn't needed, but doesn't hurt either)
		flice1_vals.push_back(std::make_tuple(lindex[0], lindex[1], ii->second));
	}

	// Zero out the GLINT2-only version completely
	auto flice1_im(flice1_im_f.to_blitz());
	flice1_im = 0;

	// Replace with our values
	for (auto ii=flice1_vals.begin(); ii != flice1_vals.end(); ++ii) {
		int ix_i = std::get<0>(*ii);
		int ix_j = std::get<1>(*ii);
		double val = std::get<2>(*ii);

		flice1_im(ix_i, ix_j) += val;
	}
	// -----------------------------------------------------
	// Balance flice against other landcover types
//	auto fgrnd1(fgrnd1_f.to_blitz());
//	auto focean1(focean1_f.to_blitz());
//	auto flake1(flake1_f.to_blitz());
// This correction is taken care of in ModelE (for now)
// See: FLUXES.f alloc_fluxes()
//	fgrnd1 = 1.0 - focean1 - flake1 - flice1;

#if 0
NcFile nc("flice1_1.nc", NcFile::Replace);
auto flice1_c(giss::f_to_c(flice1));
auto flice1_im_c(giss::f_to_c(flice1_im));
auto a(giss::netcdf_define(nc, "flice1", flice1_c));
auto b(giss::netcdf_define(nc, "flice1_im", flice1_im_c));
a();
b();
nc.close();
#endif

	// -----------------------------------------------------
printf("END glint2_modele_get_flice_im_c()\n");
}
// -----------------------------------------------------
/** Produces the (dense) FHC_IM array from the (sparse) hp_to_atm
coming from rawl Glint2. */
extern "C"
void glint2_modele_get_fhc_im_c(glint2::modele::glint2_modele *api,
	giss::F90Array<double, 3> &fhc_im1h_f)	// OUT
{
	auto fhc_im1h(fhc_im1h_f.to_blitz());

	HCIndex &hc_index(*api->gcm_coupler.maker->hc_index);
	std::unique_ptr<giss::VectorSparseMatrix> hp_to_atm(api->gcm_coupler.maker->hp_to_atm());
	ModelEDomain &domain(*api->domain);

	// Filter this array, and convert to fhc format
	for (auto ii = hp_to_atm->begin(); ii != hp_to_atm->end(); ++ii) {
		int lindex1a[domain.num_local_indices];

		// Input: HP space
		int lindex[domain.num_local_indices];
		int hp1b, i1b;
		int i3b = ii.col();
		hc_index.index_to_ik(i3b, i1b, hp1b);
		domain.global_to_local(i1b, lindex);
		if (!domain.in_domain(lindex)) {
			//printf("Not in domain: i3b=%d (%d, %d, %d)\n", i3b, lindex[0], lindex[1], hp1b);
			continue;
		}

		// Output: GCM grid
		int i1a = ii.row();
		if (i1a != i1b) {
			fprintf(stderr, "HP2ATM matrix is non-local!\n");
			giss::exit(1);
		}

		// Now fill in FHC_IM
		// +1 for C-to-Fortran conversion
		fhc_im1h(lindex[0], lindex[1], hp1b+1) += ii.val();
	}
	hp_to_atm.release();


	// In a perfect world, summing FHC over elevation points will
	// sum to one.  But in reality, it sums to something else, depending
	// on size difference between grids on sphere vs. on the plane.

}
// -----------------------------------------------------
extern "C"
void glint2_modele_get_elevhp_im_c(glint2::modele::glint2_modele *api,
	giss::F90Array<double, 3> &elev1h_f)	// IN/OUT
{

	// =================== elev1h
	// Just copy out of hpdefs array, elevation points are the same
	// on all grid cells.

	auto elev1h(elev1h_f.to_blitz());
	int nhp_im = api->gcm_coupler.maker->nhp(-1);
	if (nhp_im != elev1h.extent(2)) {
		fprintf(stderr, "glint2_modele_get_elev1h: Inconsistent nhp (%d vs %d)\n", elev1h.extent(2), nhp_im);
		giss::exit(1);
	}

	// Copy 1-D height point definitions to elev1h
	for (int k=0; k < nhp_im; ++k) {
		double val = api->gcm_coupler.maker->hpdefs[k];
		for (int j=elev1h.lbound(1); j <= elev1h.ubound(1); ++j) {
		for (int i=elev1h.lbound(0); i <= elev1h.ubound(0); ++i) {
			// +1 for C-to-Fortran conversion
			elev1h(i,j,k+1) = val;
		}}
	}
}
// -----------------------------------------------------
extern "C"
void glint2_modele_init_hp_to_ices(glint2::modele::glint2_modele *api)
{
	int const rank = api->gcm_coupler.rank();	// MPI rank; debugging
printf("[%d] BEGIN glint2_modele_init_hp_to_ices\n", rank);

	ModelEDomain &domain(*api->domain);
	HCIndex &hc_index(*api->gcm_coupler.maker->hc_index);

	// ====================== hp_to_ices
	api->hp_to_ices.clear();
	for (auto sheet=api->gcm_coupler.maker->sheets.begin(); sheet != api->gcm_coupler.maker->sheets.end(); ++sheet) {

		// Get matrix for HP2ICE
		std::unique_ptr<giss::VectorSparseMatrix> imat(
			sheet->hp_to_iceinterp(IceInterp::ICE));
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

printf("[%d] END glint2_modele_init_hp_to_ices\n", rank);
}
// -----------------------------------------------------
static void densify_gcm_inputs_onroot(glint2_modele *api,
	std::vector<giss::VectorSparseVector<int,double>> const &gcm_ivals_global,
	giss::F90Array<double,3> &gcm_inputs_d_f)
{
	// This should only be run on the MPI root node
	if (!api->gcm_coupler.am_i_root()) return;

	printf("BEGIN densify_gcm_inputs_onroot()\n");

	int const rank = api->gcm_coupler.rank();	// MPI rank; debugging
	giss::CouplingContract const &contract(api->gcm_coupler.gcm_inputs);
	int n1 = api->gcm_coupler.maker->n1();

	// Fortran-style array: i,j,ihp, indexing starts at 1
	blitz::Array<double,3> gcm_inputs_d(gcm_inputs_d_f.to_blitz());

	// We ARE the root note --- densify the data into the global gcm_inputs array
	for (long ix = 0; ix < contract.size_nounit(); ++ix) {
		int ihp = api->gcm_inputs_ihp[ix];				// First elevation point for this variable
		int var_nhp = api->gcm_inputs_ihp[ix+1] - ihp;	// # elevation points for this variable.

		// Check bounds
		if (ihp+var_nhp > gcm_inputs_d.extent(2)) {
			fprintf(stderr, "[%d] gcm_inputs_d[nhp=%d] is too small (needs at least %d)\n", rank, gcm_inputs_d.extent(2), api->gcm_inputs_ihp[contract.size_nounit()]); //ihp+var_nhp);
			giss::exit(1);
		}

		// Ignore elevation point = 0 (for ELEVATION grid only),
		// which is reserved for ModelE's "legacy" elevation point.
		int modele_ihp = ihp;
		int modele_var_nhp = var_nhp;
		if (modele_var_nhp > 1) {
			// We have an ELEVATION grid destination (not ATMOSPHERE)
			modele_ihp += 1;
			modele_var_nhp -= 1;
		}

		// Index into our big array-of-array of all gcm_inputs
		blitz::Array<double,1> dense1d(
			&gcm_inputs_d(1,1,modele_ihp),	// i,j,ihp
			blitz::shape(n1*modele_var_nhp), blitz::neverDeleteData);

		// Convert this sparse vector...
		printf("Setting gcm_input %s to %g\n", contract.field(ix).name.c_str(), contract.field(ix).default_value);
//		dense1d = contract.field(ix).default_value;
		dense1d = 0;
		for (auto ii=gcm_ivals_global[ix].begin(); ii != gcm_ivals_global[ix].end(); ++ii) {
			int const i1 = ii->first;
			double const val = ii->second;
			dense1d(i1) += val;
		}
	}

	printf("END densify_gcm_inputs_onroot()\n");


}
// -------------------------------------------------------------
/** @param hpvals Values on height-points GCM grid for various fields
	the GCM has decided to provide.

	@param gcm_inputs_d_f Global (gathered) array of the Glint2
	outputs to be fed into the GCM.  This only needs to be allocated
	if api->gcm_coupler.am_i_root().

	Inputs variables implied by repeated previous calls to glint2_modele_set_gcm_output().
*/
extern "C"
void  glint2_modele_couple_to_ice_c(
glint2_modele *api,
int itime,
giss::F90Array<double,3> &gcm_inputs_d_f)
{
	int rank = api->gcm_coupler.rank();	// MPI rank; debugging
	double time_s = itime * api->dtsrc;

	printf("[%d] BEGIN glint2_modele_couple_to_ice_c(itime=%d, time_s=%f, dtsrc=%f)\n", rank, itime, time_s, api->dtsrc);

//	GCMCoupler &coupler(api->gcm_coupler);

	if (api->gcm_coupler.gcm_out_file.length() > 0) {
		// Write out to DESM file
		save_gcm_outputs(api, time_s, api->gcm_outputs);
	}

	// Count total number of elements in the matrices
	// (_l = local to this MPI node)
	int nele_l = 0; //api->gcm_coupler.maker->ice_matrices_size();
printf("glint2_modele_couple_to_ice_c(): hp_to_ices.size() %ld\n", api->hp_to_ices.size());
	for (auto ii = api->hp_to_ices.begin(); ii != api->hp_to_ices.end(); ++ii) {
		nele_l += ii->second.size();
	}

	// Find the max. number of fields (for input) used for any ice sheet.
	// This will determine the size of our MPI messages.
	int nfields_max = 0;
	for (auto sheet=api->gcm_coupler.maker->sheets.begin(); sheet != api->gcm_coupler.maker->sheets.end(); ++sheet) {
		int sheetno = sheet->index;
		giss::VarTransformer &vt(api->gcm_coupler.models[sheetno]->var_transformer[IceModel::INPUT]);
		nfields_max = std::max(nfields_max, (int)vt.dimension(giss::VarTransformer::OUTPUTS).size_nounit());
	}

	// Allocate buffer for that amount of stuff
printf("glint2_modele_couple_to_ice_c(): nfields_max=%d, nele_l = %d\n", nfields_max, nele_l);
	giss::DynArray<SMBMsg> sbuf(SMBMsg::size(nfields_max), nele_l);

	// Fill it in by doing a sparse multiply...
	// (while translating indices to local coordinates)
	HCIndex &hc_index(*api->gcm_coupler.maker->hc_index);
	int nmsg = 0;
printf("[%d] hp_to_ices.size() = %ld\n", rank, api->hp_to_ices.size());
	for (auto ii = api->hp_to_ices.begin(); ii != api->hp_to_ices.end(); ++ii) {
		int sheetno = ii->first;
		giss::VarTransformer &vt(api->gcm_coupler.models[sheetno]->var_transformer[IceModel::INPUT]);


		std::vector<hp_to_ice_rec> &mat(ii->second);

printf("[%d] mat[sheetno=%d].size() == %ld\n", rank, sheetno, mat.size());
		// Skip if we have nothing to do for this ice sheet
		if (mat.size() == 0) continue;

		// Get the CSR sparse matrix to convert GCM outputs to ice model inputs
		giss::CSRAndUnits trans = vt.apply_scalars({
			std::make_pair("by_dt", 1.0 / ((itime - api->itime_last) * api->dtsrc)),
			std::make_pair("unit", 1.0)});

		// Do the multiplication
		for (int j=0; j < mat.size(); ++j) {
			hp_to_ice_rec &jj(mat[j]);
			SMBMsg &msg = sbuf[nmsg];
			msg.sheetno = sheetno;
			msg.i2 = jj.row;

#if 1
			// Convert from (GCM output) to (Ice Model input) units while regridding
			for (int xi=0; xi<vt.dimension(giss::VarTransformer::OUTPUTS).size_nounit(); ++xi) {
				double inp = 0;
				std::vector<std::pair<int, double>> const &row(trans.mat[xi]);
				for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
					int xj = xjj->first;
					double io_val = xjj->second;
					inp += io_val * (*api->gcm_outputs[xj])(jj.col_i, jj.col_j, jj.col_k);
				}
				msg[xi] = jj.val * (inp + trans.units[xi]);
			}
#else
			// This is the old code for the above (that doesn't convert units)
			msg[0] = jj.val * (*api->gcm_outputs[0])(jj.col_i, jj.col_j, jj.col_k);
			msg[1] = jj.val * (*api->gcm_outputs[1])(jj.col_i, jj.col_j, jj.col_k);
			msg[2] = jj.val * (*api->gcm_outputs[2])(jj.col_i, jj.col_j, jj.col_k);
#endif
//			// This is even older code (variables are hardwired)
//			msg[0] = jj.val * smb1h(jj.col_i, jj.col_j, jj.col_k);
//			msg[1] = jj.val * seb1h(jj.col_i, jj.col_j, jj.col_k);
//			msg[2] = jj.val * tg21h(jj.col_i, jj.col_j, jj.col_k);

//printf("msg = %d (i,j, hc)=(%d %d %d) i2=%d %g %g (%g %g)\n", msg.sheetno, lindex[0], lindex[1], ihp+1, msg.i2, msg[0], msg[1], smb1h(lindex[0], lindex[1], ihp+1), seb1h(lindex[0], lindex[1], ihp+1));

			++nmsg;
		}
	}

	// Sanity check: make sure we haven't overrun our buffer
	if (nmsg != sbuf.size) {
		fprintf(stderr, "Wrong number of items in buffer: %d vs %d expected\n", nmsg, sbuf.size);
		giss::exit(1);
	}

	// sbuf has elements for ALL ice sheets here
	giss::CouplingContract const &contract(api->gcm_coupler.gcm_inputs);
	std::vector<giss::VectorSparseVector<int,double>> gcm_ivals_global(contract.size_nounit());
	api->gcm_coupler.couple_to_ice(time_s, nfields_max, sbuf, gcm_ivals_global);

	// Decode the outputs to a dense array for ModelE
	int nhp = glint2_modele_nhp_gcm(api);	// == Glint2's nhp + 1, for legacy ice in ModelE
	int n1 = api->gcm_coupler.maker->n1();


	if (api->gcm_coupler.am_i_root()) {
		// auto gcm_inputs_d(gcm_inputs_d_f.to_blitz());

		densify_gcm_inputs_onroot(api, gcm_ivals_global, gcm_inputs_d_f);

		if (api->gcm_coupler.gcm_in_file.length() > 0) {
			// Write out to DESM file
			save_gcm_inputs(api, time_s, gcm_inputs_d_f);
		}
	}

	api->itime_last = itime;

	printf("[%d] END glint2_modele_couple_to_ice_c(itime=%d)\n", rank, itime);
}
// -------------------------------------------------------------
extern "C"
void  glint2_modele_get_initial_state_c(
glint2_modele *api,
giss::F90Array<double,3> &gcm_inputs_d_f)
{
	int rank = api->gcm_coupler.rank();	// MPI rank; debugging

//	auto gcm_inputs_d(gcm_inputs_d_f.to_blitz());

	printf("[%d] BEGIN glint2_modele_get_initial_state_c()\n", rank);

	GCMCoupler &coupler(api->gcm_coupler);

	// sbuf has elements for ALL ice sheets here
	giss::CouplingContract const &contract(api->gcm_coupler.gcm_inputs);
	std::vector<giss::VectorSparseVector<int,double>> gcm_ivals_global(contract.size_nounit());
	coupler.get_initial_state(gcm_ivals_global);

	// Decode the outputs to a dense array for ModelE
	if (api->gcm_coupler.am_i_root()) {

		densify_gcm_inputs_onroot(api, gcm_ivals_global, gcm_inputs_d_f);
		const double time_s = coupler.gcm_params.time_start_s;
		if (coupler.gcm_in_file.length() > 0) {
			// Write out to DESM file
			save_gcm_inputs(api, time_s, gcm_inputs_d_f);
		}
	}

	printf("[%d] END glint2_modele_get_initial_state_c()\n", rank);
}

extern"C"
int glint2_modele_set_gcm_output_c(
glint2_modele *api,
char const *field_name_f, int field_name_len,
giss::F90Array<double, 3> &arr_f)
{
	std::string const field_name(field_name_f, field_name_len);
	int field_ix = api->gcm_coupler.gcm_outputs.index(field_name);
	if (api->gcm_outputs.size() <= field_ix) api->gcm_outputs.resize(field_ix+1);
	blitz::Array<double,3> *arrp = new blitz::Array<double,3>(arr_f.to_blitz());
	api->gcm_outputs[field_ix].reset(arrp);
}

// ===============================================================
