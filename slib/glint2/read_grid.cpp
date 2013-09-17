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

#include <memory>

#include <giss/ncutil.hpp>
#include <glint2/Grid.hpp>
#include <glint2/Grid_XY.hpp>
#include <glint2/Grid_LonLat.hpp>
#include <glint2/ExchangeGrid.hpp>

namespace glint2 {

/** @param fname Name of file to load from (eg, an overlap matrix file)
@param vname Eg: "grid1" or "grid2" */
std::unique_ptr<Grid> read_grid(NcFile &nc, std::string const &vname)
{
	auto info_var = nc.get_var((vname + ".info").c_str());

	Grid::Type type = giss::parse_enum<Grid::Type>(
		giss::get_att(info_var, "type")->as_string(0));

printf("read_grid(type=%s)\n", type.str());
	std::unique_ptr<Grid> grid;
	switch(type.index()) {
		case Grid::Type::XY :
			grid.reset(new Grid_XY());
			break;
		case Grid::Type::LONLAT :
			grid.reset(new Grid_LonLat());
			break;
		case Grid::Type::EXCHANGE :
			grid.reset(new ExchangeGrid());
			break;
		default :
			grid.reset(new Grid(type));
			break;
	}

	grid->read_from_netcdf(nc, vname);
	return grid;
}

}
