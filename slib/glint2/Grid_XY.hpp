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

#include <boost/function.hpp>
#include <giss/blitz.hpp>
#include <glint2/Grid.hpp>

namespace glint2 {


/** Represents a Cartesian grid with non-equally-spaced grid cell boundaries. */
class Grid_XY : public Grid {

public:

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

	// -----------------------------------------
	Grid_XY() : Grid(Grid::Type::XY) {
		coordinates = Grid::Coordinates::XY;
		parameterization = Grid::Parameterization::L0;
	}

	long ij_to_index(int i, int j) const
		{ return j * nx() + i; }

	void index_to_ij(long index, int &i, int &j) const
	{
		int n = nx();
		j = index / n;
		i = index - j*n;
	}

	template<class T>
	blitz::Array<T,2> reshape_1d_to_xy(blitz::Array<T,1> &arr1d);

	template<class T>
	blitz::Array<T,1> reshape_xy_to_1d(blitz::Array<T,2> &arrxy);


	template<class T>
	blitz::Array<T,2> const reshape_1d_to_xy(blitz::Array<T,1> const &arr1d);

	template<class T>
	blitz::Array<T,1> const reshape_xy_to_1d(blitz::Array<T,2> const &arr_xy);



	/** Create a new Cartesian grid with arbitrary grid cell boundaries.
	@param name Value of <generic-name>.info:name in netCDF file.
	@param xb x-axis boundaries of grid cells, sorted low to high.
	@param yb y-axis boundaries of grid cells, sorted low to high.
	@param euclidian_clip Only realize grid cells that pass this test.
	@return The newly created Grid_XY.
	@see EuclidianClip
	*/
	void realize(boost::function<bool(Cell const &)> const &euclidian_clip);

	// ------------------------------------------------

//private:
//	void netcdf_write(
//		boost::function<void()> const &parent,
//		NcFile *nc, std::string const &vname);

public:
	virtual boost::function<void()> netcdf_define(NcFile &nc, std::string const &vname) const;

	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);

};

// -----------------------------------------------------------------

template<class T>
blitz::Array<T,2> Grid_XY::reshape_1d_to_xy(blitz::Array<T,1> &arr_1d)
	{ return giss::reshape<T, 1, 2>(arr_1d, blitz::shape(ny(), nx())); }

template<class T>
blitz::Array<T,1> Grid_XY::reshape_xy_to_1d(blitz::Array<T,2> &arr_xy)
	{ return giss::reshape<T, 2, 1>(arr_xy, blitz::shape(ny() * nx())); }


/// Const versions
template<class T>
blitz::Array<T,2> const Grid_XY::reshape_1d_to_xy(blitz::Array<T,1> const &arr_1d)
	{ return giss::reshape<T, 1, 2>(arr_1d, blitz::shape(ny(), nx())); }

template<class T>
blitz::Array<T,1> const Grid_XY::reshape_xy_to_1d(blitz::Array<T,2> const &arr_xy)
	{ return giss::reshape<T, 2, 1>(arr_xy, blitz::shape(ny() * nx())); }


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
@return The newly created Grid_XY.
@see EuclidianClip
*/
void set_xy_boundaries(Grid_XY &grid,
	double x0, double x1, double dx,
	double y0, double y1, double dy);

void set_xy_centers(Grid_XY &grid,
	double x0, double x1, double dx,
	double y0, double y1, double dy);



}
