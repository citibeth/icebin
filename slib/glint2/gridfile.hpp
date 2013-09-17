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

#include <string>
#include "Grid.hpp"

/** Top-level function to read a grid from a file */
std::unique_ptr<Grid> read_gridfile(std::string fname)
{
	NcFile nc(fname.c_str(), NcFile::ReadOnly);
	auto ret = Grid::netcdf_read(nc, "grid");
	nc.close();
	return ret;
}

/** Top-level function to write a grid to a file */
void write_gridfile(Grid const &grid, std::string fname)
{
	NcFile nc(fname.c_str(), NcFile::Replace);
	writefn = grid.netcdf_define(nc, "grid");
	writefn();
	nc.close();
}

