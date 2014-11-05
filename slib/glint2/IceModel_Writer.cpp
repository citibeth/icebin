/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013, 2014 by Robert Fischer
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

#include <mpi.h>		// Must be first
#include <boost/filesystem.hpp>
#include <glint2/IceModel_Writer.hpp>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <glint2/GCMParams.hpp>
#include <glint2/GCMCoupler.hpp>

namespace glint2 {


std::vector<NcDim const *> IceModel_Writer::add_dims(NcFile &nc)
{
	printf("BEGIN add_dims()\n");

	std::vector<NcDim const *> ret;
	unsigned int ndim = dim_names.size();
	ret.reserve(ndim);
	ret.push_back(nc.add_dim(dim_names[0].c_str()));	// No dimsize --> unlimited
	for (unsigned int i=1; i < ndim; ++i) {
		ret.push_back(nc.add_dim(dim_names[i].c_str(), counts[i]));
	}

	printf("END add_dims()\n");
	return ret;
}

std::vector<NcDim const *> IceModel_Writer::get_dims(NcFile &nc)
{
	std::vector<NcDim const *> ret;
	unsigned int ndim = dim_names.size();
	ret.reserve(ndim);
	for (unsigned int i=0; i < ndim; ++i)
		ret.push_back(nc.get_dim(dim_names[i].c_str()));
	return ret;
}


/** Specialized init signature for IceModel_Writer */
void IceModel_Writer::init(
	std::shared_ptr<glint2::Grid> const &grid2,
	IceModel const *model, std::string const &sheet_name)

{
	printf("BEGIN IceModel_Writer::init(%s)\n", sheet_name.c_str());

	IceModel::init(grid2);

	main_model = model;

	// Try to be clever about making multi-dimensional arrays
	// in the output according to the grid the user expects.
	switch(grid2->type.index()) {
		case Grid::Type::XY : {
			Grid_XY const *grid2_xy = dynamic_cast<Grid_XY const *>(&*grid2);
			// First dimension will be unlimited
			dim_names = {"time", "ny", "nx"};
			counts = {1, grid2_xy->ny(), grid2_xy->nx()};
			cur = {0, 0, 0};
		} break;
		default : {
			// First dimension will be unlimited
			dim_names = {"time", "n2"};
			counts = {1, grid2->ndata()};
			cur = {0, 0};
		} break;
	};

	// Only need to run one copy of this
	GCMParams const &gcm_params(coupler->gcm_params);
	if (gcm_params.gcm_rank != gcm_params.gcm_root) return;

	// Put our output files in this directory, one named per ice sheet.
	auto output_dir = boost::filesystem::absolute(
		boost::filesystem::path("ice_model_inputs"),
		coupler->gcm_params.config_dir);
	boost::filesystem::create_directory(output_dir);	// Make sure it exists

	// Set up the output file
	// Create netCDF variables based on details of the coupling contract.xs
	output_fname = (output_dir / (sheet_name + ".nc")).string();
	printf("IceModel_Writer opening file %s\n", output_fname.c_str());

	printf("END IceModel_Writer::init_from_ice_model(%s)\n", sheet_name.c_str());
}

void IceModel_Writer::start_time_set()
{
	// We just need the input (or output) coupling contract
	contract[io] = main_model->contract[io];

}

void IceModel_Writer::init_output_file()
{
	GCMParams const &gcm_params(coupler->gcm_params);

	NcFile nc(output_fname.c_str(), NcFile::Replace);
	std::vector<const NcDim *> dims = add_dims(nc);
	NcDim *one_dim = nc.add_dim("one", 1);

	NcVar *time0_var = nc.add_var("time0", giss::get_nc_type<double>(), one_dim);
	time0_var->add_att("units", gcm_params.time_units.c_str());
	time0_var->add_att("calendar", "365_day");
	time0_var->add_att("axis", "T");
	time0_var->add_att("long_name", "Simulation start time");

	NcVar *time_var = nc.add_var("time", giss::get_nc_type<double>(), dims[0]);
	time_var->add_att("units", gcm_params.time_units.c_str());
	time_var->add_att("calendar", "365_day");
	time_var->add_att("axis", "T");
	time_var->add_att("long_name", "Coupling times");


	for (auto cf = contract[io].begin(); cf != contract[io].end(); ++cf) {
		NcVar *var = nc.add_var(cf->name.c_str(), giss::get_nc_type<double>(),
			dims.size(), &dims[0]);
		var->add_att("units", cf->units.c_str());
		var->add_att("description", cf->description.c_str());
	}

	// Put initial time in it...
	long cur_b[1]{0};
	long counts_b[1]{1};
	time0_var->set_cur(cur_b);
	time0_var->put(&gcm_params.time_start_s, counts_b);

	nc.close();

	output_file_initialized = true;
}

/** @param index Index of each grid value.
@param vals The values themselves -- could be SMB, Energy, something else...
TODO: More params need to be added.  Time, return values, etc. */
void IceModel_Writer::run_decoded(double time_s,
	std::vector<blitz::Array<double,1>> const &ivals2,
	std::vector<blitz::Array<double,1>> &ovals2)	// Not used for IceModel_Writer
{
	// Only need to run one copy of this
	GCMParams const &gcm_params(coupler->gcm_params);
	if (gcm_params.gcm_rank != gcm_params.gcm_root) return;

printf("BEGIN IceModel_Writer::run_decoded\n");
	if (!output_file_initialized) init_output_file();
	NcFile nc(output_fname.c_str(), NcFile::Write);

	// Read index info
	NcDim *time_dim = nc.get_dim("time");		// No dimsize --> unlimited
	cur[0] = time_dim->size();

	// Write the current time
	NcVar *time_var = nc.get_var("time");
	time_var->set_cur(&cur[0]);
	time_var->put(&time_s, &counts[0]);


	// Write the other variables
	int i = 0;
	for (auto cf = contract[io].begin(); cf != contract[io].end(); ++cf, ++i) {
		NcVar *ncvar = nc.get_var(cf->name.c_str());
		ncvar->set_cur(&cur[0]);

		double const *data = (io == IceModel::INPUT ? ivals2[i].data() : ovals2[i].data());
		ncvar->put(data, &counts[0]);
	}

	nc.close();
printf("END IceModel_Writer::run_decoded\n");
	
}


}
