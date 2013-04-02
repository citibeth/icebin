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

