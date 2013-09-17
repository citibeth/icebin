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

#include <glint2/IceModel_DISMAL.hpp>
#include <cstdio>
#include <cmath>

namespace glint2 {

/** Query all the ice models to figure out what fields they need */
void IceModel_DISMAL::get_required_fields(std::set<IceField> &fields)
{
	fields.insert(IceField::MASS_FLUX);
	fields.insert(IceField::ENERGY_FLUX);
	fields.insert(IceField::SURFACE_T);
	fields.insert(IceField::TG2);
}



/** @param index Index of each grid value.
@param vals The values themselves -- could be SMB, Energy, something else...
TODO: More params need to be added.  Time, return values, etc. */
void IceModel_DISMAL::run_decoded(long itime,
	std::map<IceField, blitz::Array<double,1>> const &vals2)
{
	auto mass(get_field(vals2, IceField::MASS_FLUX));
	auto energy(get_field(vals2, IceField::ENERGY_FLUX));
	auto surfacet(get_field(vals2, IceField::SURFACE_T));
	auto tg2(get_field(vals2, IceField::TG2));
#if 0
	// Re-shape the arrays
	blitz::Array<double,2> const mass(
		const_cast<double *>(vals2.find(IceField::MASS_FLUX)->second.data()),
		blitz::shape(ny,nx), blitz::neverDeleteData);
	blitz::Array<double,2> const energy(
		const_cast<double *>(vals2.find(IceField::ENERGY_FLUX)->second.data()),
		blitz::shape(ny,nx), blitz::neverDeleteData);
	blitz::Array<double,2> surfacet(
		const_cast<double *>(vals2.find(IceField::SURFACE_T)->second.data()),
		blitz::shape(ny,nx), blitz::neverDeleteData);
#endif

	char fname[100];
	sprintf(fname, "dismal-%ld.nc", itime);
	NcFile ncout(fname, NcFile::Replace);

	std::vector<boost::function<void ()>> fns;
	NcDim *nx_dim = ncout.add_dim("nx", nx);
	NcDim *ny_dim = ncout.add_dim("ny", ny);

	// Define variables
	fns.push_back(giss::netcdf_define(ncout, "mass", mass, {ny_dim, nx_dim}));
	fns.push_back(giss::netcdf_define(ncout, "energy", energy, {ny_dim, nx_dim}));
	fns.push_back(giss::netcdf_define(ncout, "T", surfacet, {ny_dim, nx_dim}));
	fns.push_back(giss::netcdf_define(ncout, "TG2", tg2, {ny_dim, nx_dim}));

	// Write data to netCDF file
	for (auto ii = fns.begin(); ii != fns.end(); ++ii) (*ii)();
	ncout.close();
}


}
