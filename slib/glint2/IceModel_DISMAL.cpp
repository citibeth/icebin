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
#include <glint2/IceModel_DISMAL.hpp>
#include <cstdio>
#include <cmath>
#include <cassert>

namespace glint2 {

// Arguments that are paths, and thus need pathname resolution
static std::set<std::string> path_args = {"output_dir"};

void IceModel_DISMAL::IceModel_DISMAL::init(
		IceModel::GCMParams const &_gcm_params,
		std::shared_ptr<glint2::Grid> const &grid2,
		NcFile &nc,
		std::string const &vname_base,
		NcVar *const_var)
{
	printf("BEGIN IceModel_DISMAL::init(%s)\n", vname_base.c_str());

	IceModel_Decode::init(_gcm_params, grid2->ndata());
	auto grid2_xy = dynamic_cast<Grid_XY const *>(&*grid2);
	nx = grid2_xy->nx();
	ny = grid2_xy->ny();
	auto dismal_var = giss::get_var_safe(nc, (vname_base + ".dismal").c_str());	// DISMAL parameters
	output_dir = boost::filesystem::absolute(
		boost::filesystem::path(giss::get_att(dismal_var, "output_dir")->as_string(0)),
		gcm_params.config_dir);
	printf("END IceModel_DISMAL::int()\n");
}


blitz::Array<double,2> const IceModel_DISMAL::reshape_xy(
	blitz::Array<double,1> const &vals2)
{
	const double *data = vals2.data();
	return blitz::Array<double,2>(const_cast<double *>(data),
		blitz::shape(ny,nx), blitz::neverDeleteData);
}


/** @param index Index of each grid value.
@param vals The values themselves -- could be SMB, Energy, something else...
TODO: More params need to be added.  Time, return values, etc. */
void IceModel_DISMAL::run_decoded(double time_s,
	std::vector<blitz::Array<double,1>> const &vals2)
{
	// Only need to run one copy of this
	if (gcm_params.gcm_rank != gcm_params.gcm_root) return;

printf("BEGIN IceModel_DISMAL::run_decoded\n");

	char fname[30];
	long time_day = (int)(time_s / 86400. + .5);
	sprintf(fname, "%ld-dismal.nc", time_day);
	auto full_fname(output_dir / fname);
	printf("IceModel_DISMAL writing to: %s\n", full_fname.c_str());
	NcFile ncout(full_fname.c_str(), NcFile::Replace);
	assert(ncout.is_valid());

	std::vector<boost::function<void ()>> fns;
	NcDim *nx_dim = ncout.add_dim("nx", nx);
	assert(nx_dim);
	NcDim *ny_dim = ncout.add_dim("ny", ny);
	assert(ny_dim);

printf("contract[INPUT].size_nounit() == %ld\n", contract[INPUT].size_nounit());
	// Define variables
	for (int i=0; i < contract[INPUT].size_nounit(); ++i) {
		CoupledField const &field(contract[INPUT].field(i));
		printf("IceModel_DISMAL: Defining variable %s\n", field.name.c_str());

		// Convert from 1D indexing to 2D
//std::cout << "strides = " << vals2[i].stride() << std::endl;
//std::cout << "extents = " << vals2[i].extent() << std::endl;
		auto val2_xy(reshape_xy(vals2[i]));
		fns.push_back(giss::netcdf_define(ncout, "", field, val2_xy, {ny_dim, nx_dim}));
	}

	// Write data to netCDF file
	for (auto ii = fns.begin(); ii != fns.end(); ++ii) (*ii)();
	ncout.close();
printf("END IceModel_DISMAL::run_decoded\n");
}


}
