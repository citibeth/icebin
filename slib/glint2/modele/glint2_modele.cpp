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

using namespace glint2;
using namespace glint2::modele;


// ================================================================
// Stuff from LANDICE_DRV.f
/**
@param vals1d Base of the ij array to scatter.
@param n Number of elements in the ij array (== im * jm)
@param out_ix Specifies which scattered variable to store in, on Fortran side (see gcm_input CouplingContract). */
extern "C" void store_li_output_ij(double *vals1d, int n, int out_ix);

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
extern "C"
void glint2_modele_init_ncfile(glint2_modele *api,
giss::CouplingContract const &contract,
std::string const &fname)
{
	GCMParams const &gcm_params(api->gcm_coupler.gcm_params);

	// Set up NetCDF file to store GCM output as we received them (modele_out.nc)
	int nhp = glint2_modele_nhp(api);
	NcFile ncout(fname.c_str(), NcFile::Replace);
	NcDim *im_dim = ncout.add_dim("im", api->domain->im);
	NcDim *jm_dim = ncout.add_dim("jm", api->domain->jm);
	NcDim *nhp_dim = ncout.add_dim("nhp", nhp);
	NcDim *one_dim = ncout.add_dim("one", 1);
	NcDim *time_dim = ncout.add_dim("time");		// No dimsize --> unlimited
	const NcDim *dims[4]{time_dim, nhp_dim, jm_dim, im_dim};

	const NcDim *dims_b[1]{one_dim};
	NcVar *time0_var = ncout.add_var("time0", giss::get_nc_type<double>(), 1, dims_b);
	time0_var->add_att("units", gcm_params.time_units.c_str());
	time0_var->add_att("calendar", "365_day");
	time0_var->add_att("axis", "T");
	time0_var->add_att("long_name", "Simulation start time");

	NcVar *time_var = ncout.add_var("time", giss::get_nc_type<double>(), 1, dims);
	time_var->add_att("units", gcm_params.time_units.c_str());
	time_var->add_att("calendar", "365_day");
	time_var->add_att("axis", "T");
	time_var->add_att("long_name", "Coupling times");

	//giss::CouplingContract &gcm_outputs(api->gcm_coupler.gcm_outputs);

	for (unsigned int i=0; i < contract.size_nounit(); ++i) {
		NcVar *nc_var = ncout.add_var(contract.name(i).c_str(),
			giss::get_nc_type<double>(), 4, dims);

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
}
// -----------------------------------------------------
void glint2_modele_save_gcm_inputs(
glint2_modele *api,
double time_s,
blitz::Array<double,3> gcm_inputs)
{
	// Get dimensions of full domain
	int nhp = glint2_modele_nhp(api);
	GCMCoupler &coupler(api->gcm_coupler);
	giss::CouplingContract const &contract(coupler.gcm_inputs);

	// Open output netCDF file
	NcFile ncout(api->gcm_coupler.gcm_in_file.c_str(), NcFile::Write);	// Read/Write
	NcDim *time_dim = ncout.get_dim("time");
	NcDim *nhp_dim = ncout.get_dim("nhp");
	NcDim *jm_dim = ncout.get_dim("jm");
	NcDim *im_dim = ncout.get_dim("im");

	long cur_ijhc[4]{time_dim->size(),0,0,0};		// time, nhp, jm, im
	long counts_ijhc[4]{1, nhp_dim->size(), jm_dim->size(), im_dim->size()};

	long cur_ij[4]{time_dim->size(),0,0};		// time, nhp, jm, im
	long counts_ij[4]{1, jm_dim->size(), im_dim->size()};

	NcVar *time_var = ncout.get_var("time");
	time_var->set_cur(cur_ijhc);	// Could be cur_ij, either way is fine
	time_var->put(&time_s, counts_ijhc);

	// Write arrays to it
	int base_index = 0;
	for (unsigned int i=0; i < contract.size_nounit(); ++i) {
		double const *array_base = &gcm_inputs(base_index,0,0);

		NcVar *nc_var = ncout.get_var(contract.name(i).c_str());

		if (contract.field(i).grid == "ATMOSPHERE") {
			nc_var->set_cur(cur_ij);
			nc_var->put(array_base, counts_ij);
			base_index += 1;
		} else if (contract.field(i).grid == "ELEVATION") {
			nc_var->set_cur(cur_ijhc);
			nc_var->put(array_base, counts_ijhc);
			base_index += nhp;
		}
	}

	ncout.close();
}
// -----------------------------------------------------
/** @param hpvals Values on height-points GCM grid for various fields
	the GCM has decided to provide. */
void  glint2_modele_save_gcm_outputs(
glint2_modele *api,
double time_s,
std::vector<blitz::Array<double,3>> &inputs)
{

	// Get dimensions of full domain
	int nhp = glint2_modele_nhp(api);
	ModelEDomain const *domain(&*api->domain);

	int const rank = api->gcm_coupler.rank();	// MPI rank; debugging

	GCMCoupler &coupler(api->gcm_coupler);
	printf("[%d] BEGIN glint2_modele_save_gcm_outputs(time_s=%f)\n", rank, time_s);

	giss::CouplingContract &gcm_outputs(coupler.gcm_outputs);

	// Count total number of elements in the inputs (for this MPI domain)
	blitz::Array<double,3> &input0(inputs[0]);
	int nele_l =
		(input0.ubound(2) - input0.lbound(2) + 1) *
		(domain->j1_f - domain->j0_f + 1) *
		(domain->i1_f - domain->i0_f + 1);

	// Find the max. number of fields (for input) used for any ice sheet.
	// This will determine the size of our MPI messages.
	int nfields = gcm_outputs.size_nounit();

	// Allocate buffer for that amount of stuff
	giss::DynArray<ModelEMsg> sbuf(ModelEMsg::size(nfields), nele_l);

	// Fill it in....
	int nmsg = 0;
	for (int k=input0.lbound(2); k<=input0.ubound(2); ++k)		// nhp
	for (int j=domain->j0_f; j <= domain->j1_f; ++j)
	for (int i=domain->i0_f; i <= domain->i1_f; ++i) {
		ModelEMsg &msg = sbuf[nmsg];
		msg.i = i;
		msg.j = j;
		msg.k = k;
		for (unsigned int l=0; l<nfields; ++l) msg[l] = inputs[l](i,j,k);
		++nmsg;
	}

	// Sanity check: make sure we haven't overrun our buffer
	if (nmsg != sbuf.size) {
		fprintf(stderr, "Wrong number of items in buffer: %d vs %d expected\n", nmsg, sbuf.size);
		throw std::exception();
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
			outputs.push_back(blitz::Array<double,3>(nhp, domain->jm, domain->im));
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

		long cur[4]{time_dim->size(),0,0,0};		// time, nhp, jm, im
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
	printf("[%d] END glint2_modele_save_gcm_outputs(time_s=%f)\n", rank, time_s);
}

// ================================================================
extern "C" glint2_modele *new_glint2_modele()
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
	printf("***** BEGIN glint2_modele_init0()\n");

	// Convert Fortran arguments
	std::string glint2_config_fname(glint2_config_fname_f, glint2_config_fname_len);
	std::string maker_vname(maker_vname_f, maker_vname_len);

	// Parse directory out of glint2_config_fname
	boost::filesystem::path glint2_config_rfname = boost::filesystem::canonical(glint2_config_fname);
printf("glint2_config_rfname = %s\n", glint2_config_rfname.c_str());
	boost::filesystem::path glint2_config_dir(glint2_config_rfname.parent_path());
std::cout << "glint2_config_dir = " << glint2_config_dir << std::endl;

	// Set up parmaeters from the GCM to the ice model
	api->gcm_coupler.gcm_params = GCMParams(
		MPI_Comm_f2c(comm_f),
		root,
		glint2_config_dir);

	if (write_constants) {
		// Store constants int NetCDF file, so we can desm without ModelE.
		giss::ConstantSet &gcm_constants(api->gcm_coupler.gcm_constants);
		GCMParams const &gcm_params(api->gcm_coupler.gcm_params);
	
		NcFile nc((gcm_params.config_dir / "modele_constants.nc").c_str(), NcFile::Replace);
		gcm_constants.netcdf_define(nc, "constants");
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
	NcFile glint2_config_nc(glint2_config_rfname.c_str(), NcFile::ReadOnly);

	// Read the coupler, along with ice model proxies
	api->gcm_coupler.read_from_netcdf(glint2_config_nc, maker_vname, std::move(mdomain));
	glint2_config_nc.close();

	// Check bounds on the IceSheets, set up any state, etc.
	// This is done AFTER setup of gcm_coupler because gcm_coupler.read_from_netcdf()
	// might change the IceSheet, in certain cases.
	// (for example, if PISM is used, elev2 and mask2 will be read from related
	// PISM input file, and the version in the GLINT2 file will be ignored)
	api->gcm_coupler.maker->realize();

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
int glint2_modele_nhp(glint2_modele const *api)
{
	int ret = api->gcm_coupler.maker->nhp(-1);	// Assume all grid cells have same # EP
	// HP/HC = 1 (Fortran) reserved for legacy "non-model" ice
    // (not part of GLINT2)
	ret += 1;
//	printf("glint2_modele_nhp() returning %d\n", ret);
	return ret;
}
// -----------------------------------------------------
/** @para var_nhc Number of elevation points for this variable.
 (equal to 1 for atmosphere variables, or nhp for elevation-grid variables)
@param return: Start of this variable in the gcm_inputs_local array (Fortran 1-based index) */
extern "C"
int glint2_modele_add_gcm_input(
glint2_modele *api,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
char const *grid_f, int grid_len,
char const *long_name_f, int long_name_len)
{
	std::string field_name(field_name_f, field_name_len);
	std::string units(units_f, units_len);
	std::string grid(grid_f, grid_len);
	std::string long_name(long_name_f, long_name_len);

	int ihp = api->gcm_inputs_ihp[api->gcm_inputs_ihp.size()-1];
	int var_nhp;
	if (grid == "ATMOSPHERE") var_nhp = 1;
	else if (grid == "ELEVATION") var_nhp = glint2_modele_nhp(api);
	else {
		fprintf(stderr, "Unrecognized grid: %s\n", grid.c_str());
		throw std::exception();
	}


	api->gcm_coupler.gcm_inputs.add_field(field_name, units, grid, long_name);

	api->gcm_inputs_ihp.push_back(ihp + var_nhp);
	int ret = ihp+1;

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

printf("BEGIN glint2_modele_set_start_time: iyear1=%d, itimei=%d, dtsrc=%f, time0_s=%f\n", iyear1, itimei, api->dtsrc, time0_s);
	api->gcm_coupler.set_start_time(
		giss::time::tm(iyear1, 1, 1),
		time0_s);

	api->itime_last = itimei;

	// Finish initialization...
	// -------------------------------------------------
	if (api->gcm_coupler.gcm_out_file.length() > 0) {
		glint2_modele_init_ncfile(api,
			api->gcm_coupler.gcm_outputs,
			api->gcm_coupler.gcm_out_file);
	}
	if (api->gcm_coupler.gcm_in_file.length() > 0) {
		glint2_modele_init_ncfile(api,
			api->gcm_coupler.gcm_inputs,
			api->gcm_coupler.gcm_in_file);
	}

printf("END glint2_modele_set_start_time\n");
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
std::cout << "fgice1_f: " << fgice1_f << std::endl;
std::cout << "fgice1_glint2_f: " << fgice1_glint2_f << std::endl;

	ModelEDomain &domain(*api->domain);

	// Reconstruct arrays, using Fortran conventions
	// (smallest stride first, whatever-based indexing it came with)

	// Get the sparse vector values
	giss::VectorSparseVector<std::pair<int,int>,double> fhc1h_s;
	giss::VectorSparseVector<int,double> fgice1_s;
	api->gcm_coupler.maker->fgice(fgice1_s);

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
// This correction is taken care of in ModelE (for now)
// See: FLUXES.f
//	fgrnd1 = 1.0 - focean1 - flake1 - fgice1;

#if 0
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
#if 0
// NOT USED

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
#endif
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
	giss::F90Array<double, 3> &elev1h_f,	// IN/OUT
	int const i0, int const j0, int const i1, int const j1)			// Array bound to write in
{
printf("init_landice_com_part2 1\n");

	// =================== elev1h
	// Just copy out of hpdefs array, elevation points are the same
	// on all grid cells.

	auto elev1h(elev1h_f.to_blitz());
	int nhp_glint2 = api->gcm_coupler.maker->nhp(-1);
	int nhp = api->gcm_coupler.maker->nhp(-1) + 1;	// Add non-model HP
	if (nhp != elev1h.extent(2)) {
		fprintf(stderr, "glint2_modele_get_elev1h: Inconsistent nhp (%d vs %d)\n", elev1h.extent(2), nhp);
		throw std::exception();
	}

	// Copy 1-D height point definitions to elev1h
	for (int k=0; k < nhp_glint2; ++k) {
		double val = api->gcm_coupler.maker->hpdefs[k];
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
		// (This allows us to easily compare running with/without height points)
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

		// Set everything from mink to maxk (inclusive) as used
		for (int k=mink; k<=maxk; ++k) used1h(i,j,k) = 1;
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

printf("END glint2_modele_init_hp_to_ices\n");
}
// -----------------------------------------------------
/** @param hpvals Values on height-points GCM grid for various fields
	the GCM has decided to provide.
	@param gcm_inputs_d_f Global (gathered) array of the Glint2
	outputs to be fed into the GCM. */

extern "C"
void  glint2_modele_couple_to_ice_c(
glint2_modele *api,
int itime,
giss::F90Array<double,3> &smb1h_f,		// kg/m^2
giss::F90Array<double,3> &seb1h_f,		// J/m^2: Latent Heat
giss::F90Array<double,3> &tg21h_f,		// C
giss::F90Array<double,3> &gcm_inputs_d_f)
{
	int rank = api->gcm_coupler.rank();	// MPI rank; debugging
	double time_s = itime * api->dtsrc;

	printf("BEGIN glint2_modele_couple_to_ice_c(itime=%d, time_s=%f, dtsrc=%f)\n", itime, time_s, api->dtsrc);

	GCMCoupler &coupler(api->gcm_coupler);
	giss::CouplingContract &gcm_outputs_contract(coupler.gcm_outputs);

	// Construct vector of GCM input arrays --- to be converted to inputs for GLINT2
	std::vector<blitz::Array<double,3>> inputs(gcm_outputs_contract.size_nounit());
	inputs[gcm_outputs_contract.index("lismb")].reference(smb1h_f.to_blitz());
	inputs[gcm_outputs_contract.index("liseb")].reference(seb1h_f.to_blitz());
	inputs[gcm_outputs_contract.index("litg2")].reference(tg21h_f.to_blitz());

	if (coupler.gcm_out_file.length() > 0) {
		// Write out to DESM file
		glint2_modele_save_gcm_outputs(api, time_s, inputs);
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
			// Convert from input to output units while regridding
			for (int xi=0; xi<vt.dimension(giss::VarTransformer::OUTPUTS).size_nounit(); ++xi) {
				double inp = 0;
				std::vector<std::pair<int, double>> const &row(trans.mat[xi]);
				for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
					int xj = xjj->first;
					double io_val = xjj->second;
					inp += io_val * inputs[xj](jj.col_i, jj.col_j, jj.col_k);
				}
				msg[xi] = jj.val * (inp + trans.units[xi]);
			}
#else
			// This is the old code for the above (that doesn't convert units)
			msg[0] = jj.val * inputs[0](jj.col_i, jj.col_j, jj.col_k);
			msg[1] = jj.val * inputs[1](jj.col_i, jj.col_j, jj.col_k);
			msg[2] = jj.val * inputs[2](jj.col_i, jj.col_j, jj.col_k);
#endif
//			// This is even older code (variables are hardwired)
//			msg[0] = jj.val * smb1h(jj.col_i, jj.col_j, jj.col_k);
//			msg[1] = jj.val * seb1h(jj.col_i, jj.col_j, jj.col_k);
//			msg[2] = jj.val * tg21h(jj.col_i, jj.col_j, jj.col_k);

//printf("msg = %d (i,j, hc)=(%d %d %d) i2=%d %g %g (%g %g)\n", msg.sheetno, lindex[0], lindex[1], ihc+1, msg.i2, msg[0], msg[1], smb1h(lindex[0], lindex[1], ihc+1), seb1h(lindex[0], lindex[1], ihc+1));

			++nmsg;
		}
	}

	// Sanity check: make sure we haven't overrun our buffer
	if (nmsg != sbuf.size) {
		fprintf(stderr, "Wrong number of items in buffer: %d vs %d expected\n", nmsg, sbuf.size);
		throw std::exception();
	}

printf("glint2_modele_couple_to_ice_c(): itime=%d, time_s=%f (dtsrc=%f)\n", itime, time_s, api->dtsrc);
	// sbuf has elements for ALL ice sheets here
	giss::CouplingContract const &contract(api->gcm_coupler.gcm_inputs);
	std::vector<giss::VectorSparseVector<int,double>> gcm_ivals_global(contract.size());
	coupler.couple_to_ice(time_s, nfields_max, sbuf, gcm_ivals_global);

	// Decode the outputs to a dense array for ModelE
	int nhp = glint2_modele_nhp(api);	// == Glint2's nhp + 1, for legacy ice in ModelE
	int n1 = api->gcm_coupler.maker->n1();


	if (api->gcm_coupler.am_i_root()) {
		blitz::Array<double,3> gcm_inputs_d(gcm_inputs_d_f.to_blitz());

		// We ARE the root note --- densify the data into the global gcm_inputs array
		for (long ix = 0; ix < contract.size(); ++ix) {
			int ihp = api->gcm_inputs_ihp[ix];
			int var_nhp = api->gcm_inputs_ihp[ix+1] - ihp;

			// Check bounds
			if (ihp+var_nhp >= gcm_inputs_d.extent(0)) {
				fprintf(stderr, "gcm_inputs_d[nhp=%d] is too small (needs at least %d)\n", gcm_inputs_d.extent(0), ihp+var_nhp);
				throw std::exception();
			}

			// Index into our big array-of-array of all gcm_inputs
			blitz::Array<double,1> dense1d(
				&gcm_inputs_d(ihp, 0, 0),
				blitz::shape(n1*var_nhp), blitz::neverDeleteData);
	
			// Convert this sparse vector...
			for (auto ii=gcm_ivals_global[ix].begin(); ii != gcm_ivals_global[ix].end(); ++ii) {
				int const i1 = ii->first;
				double const val = ii->second;
				dense1d(i1) = val;
			}
		}

		if (coupler.gcm_in_file.length() > 0) {
			// Write out to DESM file
			glint2_modele_save_gcm_inputs(api, time_s, gcm_inputs_d);
		}
	}

	api->itime_last = itime;

	printf("END glint2_modele_couple_to_ice_c(itime=%d)\n", itime);
}

// ===============================================================
