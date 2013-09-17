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

#include <memory>
#include <glint2/Grid.hpp>
#include <giss/Proj.hpp>

namespace glint2 {

class ExchangeGrid : public Grid {
public :
	/** Transformation to get from local Grid coords to Exchange Grid coords */
	giss::Proj2 proj1, proj2;

	/** Keep track of the "full" indexing space for the Overlap Matrix. */
	long grid1_ncells_full;
	long grid2_ncells_full;

	/** @param proj Projection to use to project Lon/Lat grids to XY,
	if no projection is found in the XY-type grid. */
	ExchangeGrid(Grid const &grid1, Grid const &grid2, std::string const &_sproj="");

	ExchangeGrid(): Grid(Grid::Type::EXCHANGE) {}

	virtual boost::function<void()> netcdf_define(NcFile &nc, std::string const &vname) const;

	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);

};


/** @param grid2 Put in an RTree */
extern std::unique_ptr<Grid> compute_exchange_grid(Grid &grid1, Grid &grid2);

}	// namespace glint2
