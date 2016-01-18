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
#include <functional>
#include <ibmisc/blitz.hpp>
#include <icebin/Grid.hpp>
#include <icebin/Indexing.hpp>

namespace icebin {


/** Represents a Cartesian grid with non-equally-spaced grid cell boundaries. */
struct GridSpec_XY {

	std::string name;
	std::string sproj;
	std::function<bool(Cell const &)> euclidian_clip;
	std::unique_ptr<icebin::Indexing> indexing;

	/** Cell boundaries in the x direction.
	Sorted low to high.
	Number of grid cells in the x direction = x_boundaries.size() - 1. */
	std::vector<double> xb;

	/** Cell boundaries in the y direction.
	Sorted low to high.
	Number of grid cells in the y direction = y_boundaries.size() - 1. */
	std::vector<double> yb;

	int nx() const { return xb.size() - 1; }
	int ny() const { return yb.size() - 1; }

	void make_grid(Grid &grid);
	void ncio(ibmisc::NcIO &ncio, std::string const &vname);

};

// -----------------------------------------------------------------


// -----------------------------------------------------------------

/** Create a new Cartesian grid with evenly spaced grid cell boundaries.
@param name Value of <generic-name>.info:name in netCDF file.
@param x0 Lowest boundary in the x direction.
@param x1 Highest boundary in the x direction.	
@param dx Size of grid cell in the x direction.
	Will be adjusted if (x1-x0) is not an even multiple of dx
@param y0 Lowest boundary in the y direction.
@param y1 Highest boundary in the y direction.	
@param dy Size of grid cell in the y direction.
	Will be adjusted if (y1-y0) is not an even multiple of dy
@param euclidian_clip Only realize grid cells that pass this test.
@return The newly created GridSpec_XY.
@see EuclidianClip
*/
void set_xy_boundaries(GridSpec_XY &grid,
	double x0, double x1, double dx,
	double y0, double y1, double dy);

void set_xy_centers(GridSpec_XY &grid,
	double x0, double x1, double dx,
	double y0, double y1, double dy);



}
