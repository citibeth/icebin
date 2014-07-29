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

class IceModel_Writer : public IceModel_Decode
{
	// Dimensions to use when writing to netCDF
	std::vector<std::string> dim_names;
	std::vector<long> cur;		// Base index to write in netCDF
	std::vector<long> counts;

	// The output file we are writing to...
	std::string output_fname;

public:

	IceModel_Writer(GCMCoupler const *_coupler) : IceModel_Decode(IceModel::Type::WRITER, _coupler) {}

	/** Initialize any grid information, etc. from the IceSheet struct.
	@param vname_base Construct variable name from this, out of which to pull parameters from netCDF */
	void init(
		std::shared_ptr<glint2::Grid> const &grid2,
		NcFile &nc,
		std::string const &vname_base);

	/** Init from the parent ice model.  Copy contracts, etc.  This also opens a file to write to. */
	void init_from_ice_model(IceModel const *model, std::string const &sheet_name);

	/** @param index Index of each grid value.
	@param vals The values themselves -- could be SMB, Energy, something else...
	TODO: More params need to be added.  Time, return values, etc. */
	void run_decoded(double time_s,
		std::vector<blitz::Array<double,1>> const &vals2);

protected:

	std::vector<NcDim const *> add_dims(NcFile &nc);
	std::vector<NcDim const *> get_dims(NcFile &nc);

};

}
