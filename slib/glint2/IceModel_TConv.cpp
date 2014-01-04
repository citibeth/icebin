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

#include <mpi.h>		// Must be first
#include <glint2/IceModel_TConv.hpp>

namespace glint2 {

	IceModel_TConv::IceModel_TConv(std::unique_ptr<IceModel_Decode> &&_model,
		double _LHM, double _SHI) :
		model(std::move(_model)), 
		LHM(_LHM), SHI(_SHI)
	{
		IceModel_Decode::init(_model->ndata());
	}

 	/** Query all the ice models to figure out what fields they need */
	void IceModel_TConv::get_required_fields(std::set<IceField> &fields)
	{
printf("BEGIN IceModel_TConv::get_required_fields()\n");
		model->get_required_fields(fields);
		fields.erase(IceField::SURFACE_T);
		fields.insert(IceField::MASS_FLUX);
		fields.insert(IceField::ENERGY_FLUX);
printf("END IceModel_TConv::get_required_fields()\n");
	}

	/** @param index Index of each grid value.
	@param vals The values themselves -- could be SMB, Energy, something else...
	TODO: More params need to be added.  Time, return values, etc.
	@param time_s Time since start of simulation, in seconds */
	void IceModel_TConv::run_decoded(double time_s,
		std::map<IceField, blitz::Array<double,1>> const &vals2)
	{
printf("BEGIN IceModel_TConv::run_decoded(%f)\n", time_s);

		// Copy existing fields
		std::map<IceField, blitz::Array<double,1>> ovals;
		for (auto ii = vals2.begin(); ii != vals2.end(); ++ii)
			ovals.insert(*ii);

		// Augment with SURFACE_T
		blitz::Array<double,1> mass(vals2.find(IceField::MASS_FLUX)->second);
		blitz::Array<double,1> energy(vals2.find(IceField::ENERGY_FLUX)->second);
		blitz::Array<double,1> surfacet(ndata());

		for (int i=0; i<ndata(); ++i) {
			surfacet(i) = (energy(i) / mass(i) + LHM) / SHI;
		}
		ovals.insert(std::make_pair(IceField::SURFACE_T, surfacet));

		model->run_decoded(time_s, ovals);
printf("END IceModel_TConv::run_decoded(%ld)\n", time_s);
	}

	void IceModel_TConv::read_from_netcdf(NcFile &nc, std::string const &vname)
	{
		model->read_from_netcdf(nc, vname);
	}

}
