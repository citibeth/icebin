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

#pragma once

#include <boost/enum.hpp>
#include <boost/filesystem/path.hpp>
#include <glint2/IceSheet.hpp>

namespace glint2 {

/** The different things we can pass to an ice model. */
BOOST_ENUM_VALUES( IceField, int,
	(MASS_FLUX) 	(0)		// kg/(s m^2)
	(ENERGY_FLUX)	(1)		// W/m^2
	(TG2)			(2)		// C (Mean T at bottom of firn/snow model)
	(SURFACE_T)		(3)		// C (Computed T based on mass & energy flux)
);

class IceModel {
public:
	BOOST_ENUM_VALUES( Type, int,
		(DISMAL)		(0)		// Demo Ice Sheet Model and LandIce
		(PISM)			(1)
		(ISSM)			(2)
	);

	virtual ~IceModel() {}

	/** Initialize any grid information, etc. from the IceSheet struct. */
//	virtual void init(IceSheet *sheet);

 	/** Query all the ice models to figure out what fields they need */
	virtual void get_required_fields(std::set<IceField> &fields) = 0;

	/** @param index Index of each grid value.
	@param vals The values themselves -- could be SMB, Energy, something else...
	TODO: More params need to be added.  Time, return values, etc.
	@param time_s Seconds since beginning of simulation.
	Helps with debugging. */
	virtual void run_timestep(double time_s,
		blitz::Array<int,1> const &indices,
		std::map<IceField, blitz::Array<double,1>> const &vals2) = 0;

	virtual void read_from_netcdf(NcFile &nc, std::string const &vname) {}
};

extern std::unique_ptr<IceModel> read_icemodel(
	MPI_Comm gcm_comm,
	boost::filesystem::path const &config_dir,
	NcFile &nc, std::string const &vname, IceSheet const *sheet = NULL);

}
