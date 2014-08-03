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
	/** The IceModel we're affiliated with */
	IceModel const *main_model;

	// Dimensions to use when writing to netCDF
	std::vector<std::string> dim_names;
	std::vector<long> cur;		// Base index to write in netCDF
	std::vector<long> counts;

	// The output file we are writing to...
	std::string output_fname;

public:

	IceModel_Writer(std::string const &_name, GCMCoupler const *_coupler) : IceModel_Decode(IceModel::Type::WRITER, _name, _coupler), output_file_initialized(false) {}


	/** Specialized init signature for IceModel_Writer */
	void init(
		std::shared_ptr<glint2::Grid> const &grid2,
		IceModel const *model, std::string const &sheet_name);

	void start_time_set();

protected:
	bool output_file_initialized;
	/** This is called on-demand, the first time through run_decoded(). */
	void init_output_file();

public:
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
