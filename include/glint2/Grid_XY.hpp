#pragma once

#include <boost/function.hpp>
#include "Grid.hpp"
#include "geometry.hpp"

namespace glint2 {


/** Represents a Cartesian grid with non-equally-spaced grid cell boundaries. */
class Grid_XY : public Grid {

public:

	/** Cell boundaries in the x direction.
	Sorted low to high.
	Number of grid cells in the x direction = x_boundaries.size() - 1. */
	std::vector<double> x_boundaries;

	/** Cell boundaries in the y direction.
	Sorted low to high.
	Number of grid cells in the y direction = y_boundaries.size() - 1. */
	std::vector<double> y_boundaries;


	// -----------------------------------------
	Grid_XY() : Grid("xy") {}


	void set_xy_boundaries(
		std::vector<double> const &&xb,
		std::vector<double> const &&yb)


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
	void set_xy_boundaries(
		double x0, double x1, double dx,
		double y0, double y1, double dy);

	void set_xy_centers(
		double x0, double x1, double dx,
		double y0, double y1, double dy);

	/** Create a new Cartesian grid with arbitrary grid cell boundaries.
	@param name Value of <generic-name>.info:name in netCDF file.
	@param xb x-axis boundaries of grid cells, sorted low to high.
	@param yb y-axis boundaries of grid cells, sorted low to high.
	@param euclidian_clip Only realize grid cells that pass this test.
	@return The newly created Grid_XY.
	@see EuclidianClip
	*/
	void realize_grid(boost::function<bool(Cell const &)> const &euclidian_clip);

	// ------------------------------------------------


	 boost::function<void()> netcdf_define(NcFile &nc, std::string const &generic_name) const;

	virtual std::unique_ptr<Grid::SmoothingFunction> get_smoothing_function(std::set<int> const &mask, double phi);

	void read_from_netcdf(NcFile &nc, std::string const &grid_var_name);

};

}
