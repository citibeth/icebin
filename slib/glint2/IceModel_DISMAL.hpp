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

#include <glint2/IceModel_Decode.hpp>
#include <glint2/Grid_XY.hpp>

namespace glint2 {

class IceModel_DISMAL : public IceModel_Decode
{
	// Assume a simple X-Y Cartesian grid
	int nx, ny;

public:

	IceModel_DISMAL(Grid_XY const &grid) :
		IceModel_Decode(grid),
		nx(grid.nx()), ny(grid.ny()) {}

	/** Query all the ice models to figure out what fields they need */
	void get_required_fields(std::set<IceField> &fields);

	/** @param index Index of each grid value.
	@param vals The values themselves -- could be SMB, Energy, something else...
	TODO: More params need to be added.  Time, return values, etc. */
	void run_decoded(long itime,
		std::map<IceField, blitz::Array<double,1>> const &vals2);

protected :
	blitz::Array<double,2> const get_field(
		std::map<IceField, blitz::Array<double,1>> const &vals2,
		IceField field)
	{
		const double *data = vals2.find(field)->second.data();
		return blitz::Array<double,2>(const_cast<double *>(data),
			blitz::shape(ny,nx), blitz::neverDeleteData);
	}

};

}
