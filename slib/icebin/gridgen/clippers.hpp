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

#include <icebin/gridgen/cgal.hpp>

namespace icebin {

/** A set of clipping functions to clip on the surface of a sphere.
Clipping is used to decide which grid cells in a Grid should not be realized.
Grep for "spherical_clip" to see where this might be used.
Generally, these methods need to be wrapped via boost::bind() when used.
@see Grid_LatLon */
class SphericalClip {
public:
	/** Clips everything too far from a central point.
	@param center_lon Longitude of the center point (degrees)
	@param center_lat Latitude of the center point (degrees)
	@param clip_distance_deg maximum distance from the center point for which to keep grid cells
		(degrees along the surface of a hypothetical spherical earth).
	@param lon0 Free parameter, not to be bound by bost::bind().
	@param lon1 Free parameter, not to be bound by bost::bind().
	@param lat0 Free parameter, not to be bound by bost::bind().
	@param lat1 Free parameter, not to be bound by bost::bind(). */
	static bool azimuthal(
		double center_lon, double center_lat, double clip_distance_deg,
		double lon0, double lat0, double lon1, double lat1);



	/** Clips everything outside of a lat-lon box
	@param min_lon Longitude of west edge of the clipping box (degrees).
	@param max_lon Longitude of east edge of the clipping box (degrees).
		<b>NOTE:</b> max_lon - min_lon must always be positive.  If it is not, add 360 to max_lon.
	@param min_lat Latitude of south edge of the clipping box (degrees, south pole = -90).
	@param max_lat Latitude of south edge of the clipping box (degrees, north pole = 90).
	@param lon0 Free parameter, not to be bound by bost::bind().
	@param lon1 Free parameter, not to be bound by bost::bind().
	@param lat0 Free parameter, not to be bound by bost::bind().
	@param lat1 Free parameter, not to be bound by bost::bind(). */
	static bool lonlat(
		double min_lon, double min_lat, double max_lon, double max_lat,
		double lon0, double lat0, double lon1, double lat1);

	/** Keeps all grid cells, clips nothing.
	@param lon0 Free parameter, not to be bound by bost::bind().
	@param lon1 Free parameter, not to be bound by bost::bind().
	@param lat0 Free parameter, not to be bound by bost::bind().
	@param lat1 Free parameter, not to be bound by bost::bind(). */
	static bool keep_all(double lon0, double lat0, double lon1, double lat1);

};


/** A set of clipping functions to clip in the plane.
Clipping is used to decide which grid cells in a Grid should not be realized.
Grep for "euclidian_clip" to see where this might be used.
Generally, these methods need to be wrapped via boost::bind() when used.
@see Grid_XY, Grid_LatLon */
class EuclidianClip {
public:

	/** @param clip_poly Only realize grid cells that intersect
	with this polygon (on the map)
	@param grid_cell Free parameter, not to be bound by boost::bind().
	@return True if this grid cell should be kept. */
	static bool poly(gc::Polygon_2 const &clip_poly,
		Cell const &grid_cell);

	/** Keeps all grid cells, clips nothing.
	@param grid_cell Free parameter, not to be bound by boost::bind().
	@return True */
	static bool keep_all(Cell const &grid_cell);

};

}
