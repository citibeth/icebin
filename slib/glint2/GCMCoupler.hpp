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

#include <giss/Dict.hpp>
#include <glint2/IceModel.hpp>

namespace glint2 {

class GCMCoupler {
public:

	// Only needed by root MPI node in MPI version
	giss::MapDict<int,IceModel> models;

	/** Query all the ice models to figure out what fields they need */
	std::set<IceField> get_required_fields();

	/** @param sheets (OPTIONAL): IceSheet data structures w/ grids, etc. */
	virtual void read_from_netcdf(NcFile &nc, std::string const &vname,
		std::vector<std::string> const &sheet_names,
	    giss::MapDict<std::string, IceSheet> const &sheets);

	/** Returns a unique rank number for each node in the parallel computation.
	Useful for debugging-type output. */
	virtual int rank() = 0;

};

}
