#pragma once

#include <boost/function.hpp>
#include "Proj.hpp"
#include "Grid.hpp"

namespace glint2 {

/** Represents a lat-lon grid on a sphere with non-equally-spaced grid cell boundaries.
Can optionally represent circular "cap" cells at the north and south pole,
if needed. */
class Grid_LonLat : public Grid {
	
	/** Longitude of cell boundaries (degrees), sorted low to high.
	<b>NOTE:</b> lon_boundares.last() = 360.0 + lonb.first() */
	std::vector<double> lonb;

	/** Latitude of cell boundaries (degrees), sorted low to high.
	90 = north pole, -90 = south pole. */
	std::vector<double> latb;

	/** Number of grid cell indices in longitude dimension */
	int nlon;
	/** Number of grid cell indices in latitude dimension */
	int nlat;

public:
	int nlon() { return _nlon; }
	int nlat() { return _nlat; }

	/** Number of segments used to represent each side of a grid cell with a polygonal approximation.
	Increasing this number will decrease geometric error at the
	expense of computation time for the overlap matrix. */
	int points_in_side;

	/** Projection used to make this grid.
	The projection is used to transform (spherical) lat-lon grid cells into Cartesian grid cells in the plane that can be represented with Cartesian polygons.
	<b>NOTE:</b> For best accuracy in Snowdrift, this should be an
		equal-area projection (eg: Lambert Azimuthal Equal Area). */
	Proj proj;							// Projection used to create this grid

	/** True if this grid contains a circular cap on the south pole.
	If not, then many triangular grid cells will meet at the south pole. */
	bool south_pole;

	/** True if this grid contains a circular cap on the north pole.
	If not, then many triangular grid cells will meet at the north pole. */
	bool north_pole;

	/** Radius of Earth (m), or whatever planet you're looking at.
	Used to compute theoretical exact areas of graticules. That's different
	from the radius (or elliptical shape) used in the projection. */
	double eq_rad;

	// ---------------------------------
	Grid_LonLat() : Grid("lonlat") {}

	void set_lonlatb(
		std::vector<double> const &&lonb,
		std::vector<double> const &&latb,
		bool _south_pole, bool _north_pole);


	/** Set up a Grid_LonLat with a given specification.
	@param euclidian_clip Only realize grid cells that pass this test (after projection).
	@param spherical_clip Only realize grid cells that pass this test (before projection).
	@see EuclidianClip, SphericalClip
	*/
	void realize_grid(
		boost::function<bool(double, double, double, double)> const &spherical_clip,
		boost::function<bool(Cell const &)> const &euclidian_clip);

	virtual boost::function<void()> netcdf_define(NcFile &nc, std::string const &vname) const;

	std::unique_ptr<Grid_LonLat> read_from_netcdf(NcFile &nc, std::string const &vname) const;
		
};


}
