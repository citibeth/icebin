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

#include <glint2/IceModel.hpp>

namespace glint2 {

class IceModel_Decode : public IceModel {
public:
	// Dimensionality Ice Model's vector space
	int const ndata;

	IceModel_Decode(Grid const &grid) : ndata(grid.ndata()) {}
	IceModel_Decode(int _ndata) : ndata(_ndata) {}

	/** @param index Index of each grid value.
	@param vals The values themselves -- could be SMB, Energy, something else...
	TODO: More params need to be added.  Time, return values, etc.
	@param itime Some kind of representation of the current GCM timestep.
	Helps with debugging. */
	virtual void run_timestep(long itime,
		blitz::Array<int,1> const &indices,
		std::map<IceField, blitz::Array<double,1>> const &vals2);

	/** Runs a timestep after fields have been decoded.  This is what
	one will normally want to override, unless you wish to decode
	yourself. */
	virtual void run_decoded(long itime,
		std::map<IceField, blitz::Array<double,1>> const &vals2) = 0;

};

}
