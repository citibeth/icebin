#pragma once

#include <boost/function.hpp>
#include <giss/Proj.hpp>
#include <glint2/Grid.hpp>

namespace glint2 {

/** Represents a lat-lon grid on a sphere with non-equally-spaced grid cell boundaries.
Can optionally represent circular "cap" cells at the north and south pole,
if needed. */
class Grid_LonLat : public Grid {

public:
	// ---------------------------------------------------------
	// These variables get set in concert....

	/** Longitude of cell boundaries (degrees), sorted low to high.
	<b>NOTE:</b> lon_boundares.last() = 360.0 + lonb.first() */
	std::vector<double> lonb;

	/** Latitude of cell boundaries (degrees), sorted low to high.
	90 = north pole, -90 = south pole. */
	std::vector<double> latb;

	/** True if this grid contains a circular cap on the south pole.
	If not, then many triangular grid cells will meet at the south pole. */
	bool south_pole;

	/** True if this grid contains a circular cap on the north pole.
	If not, then many triangular grid cells will meet at the north pole. */
	bool north_pole;

#if 1
	/** Number of grid cell indices in longitude dimension */
	int nlon() const { return lonb.size() - 1; }

	/** Number of grid cell indices in latitude dimension */
	int nlat() const {
		const int south_pole_offset = (south_pole ? 1 : 0);
		const int north_pole_offset = (north_pole ? 1 : 0);
		return latb.size() - 1 + south_pole_offset + north_pole_offset;
	}
#endif

	// ------------------------------------------------------

#if 0
	/** Longitude centers of grid cells */
	std::vector<double> lonc;
	std::vector<double> latc;

	/** Number of grid cell indices in longitude dimension */
	int nlon() { return lonc.size(); }

	/** Number of grid cell indices in latitude dimension */
	int nlat() { return latc.size(); }
#endif

#if 0
protected :
	int _nlon;
	int _nlat;

public :
	/** Number of grid cell indices in longitude dimension */
	int nlon() { return _nlon; }

	/** Number of grid cell indices in latitude dimension */
	int nlat() { return _nlat; }
#endif

	// ------------------------------------------------------

	/** Number of segments used to represent each side of a grid cell with a polygonal approximation.
	Increasing this number will decrease geometric error at the
	expense of computation time for the overlap matrix. */
	int points_in_side;

	/** Projection used to make this grid.
	The projection is used to transform (spherical) lat-lon grid cells into Cartesian grid cells in the plane that can be represented with Cartesian polygons.
	<b>NOTE:</b> For best accuracy in Snowdrift, this should be an
		equal-area projection (eg: Lambert Azimuthal Equal Area). */
	giss::Proj proj;							// Projection used to create this grid

	/** Radius of Earth (m), or whatever planet you're looking at.
	Used to compute theoretical exact areas of graticules. That's different
	from the radius (or elliptical shape) used in the projection. */
	double eq_rad;

	// ---------------------------------
	Grid_LonLat() : Grid("lonlat") {}

	/** This is the indexing scheme used in ModelE. */
	virtual int ij_to_index(int i, int j) const
		{ return j * nlon() + i; }

	virtual void index_to_ij(int index, int &i, int &j) const
	{
		int n = nlon();
		j = index / n;
		i = index - j*n;
	}


	/** Set up a Grid_LonLat with a given specification.
	@param euclidian_clip Only realize grid cells that pass this test (after projection).
	@param spherical_clip Only realize grid cells that pass this test (before projection).
	@see EuclidianClip, SphericalClip
	*/
	void realize(
		boost::function<bool(double, double, double, double)> const &spherical_clip,
		boost::function<bool(Cell const &)> const &euclidian_clip);

	virtual boost::function<void()> netcdf_define(NcFile &nc, std::string const &vname) const;

	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);

};


}
