#include <glint2/clippers.hpp>
#include <giss/geodesy.hpp>

namespace giss {

// We use lon0...lat1 as a poor man's substute for real CG on a spherical surface

bool SphericalClip::azimuthal(
double center_lon, double center_lat, double clip_distance_deg,
double lon0, double lat0, double lon1, double lat1)
{


//printf("distance=%g, clip_distance=%g\n", haversine_distance(center_lon, center_lat, lon0, lat0), clip_distance_deg);

	// This avoids projection-caused grid cell degeneracy
	// far away from the central area we're interested in
	if (haversine_distance(center_lon, center_lat, lon0, lat0) <= clip_distance_deg) return true;
	if (haversine_distance(center_lon, center_lat, lon1, lat0) <= clip_distance_deg) return true;
	if (haversine_distance(center_lon, center_lat, lon1, lat1) <= clip_distance_deg) return true;
	if (haversine_distance(center_lon, center_lat, lon0, lat1) <= clip_distance_deg) return true;
	return false;
}

bool SphericalClip::latlon(
double min_lon, double min_lat, double max_lon, double max_lat,
double lon0, double lat0, double lon1, double lat1)
{
	if (lon0 < min_lon && lon1 < min_lon) return false;
	if (lon0 > max_lon && lon1 > max_lon) return false;
	if (lat0 < min_lat && lat1 < min_lat) return false;
	if (lat0 > max_lat && lat1 > max_lat) return false;
	return true;
}

bool SphericalClip::keep_all(
double lon0, double lat0, double lon1, double lat1)
	{ return true; }


// =============================================================
/** @param clip_poly Only realize grid cells that intersect with this polygon (on the map) */
bool EuclidianClip::poly(gc::Polygon_2 &clip_poly,
		Cell const &grid_cell)
{
	auto cell_poly = Cell_to_Polygon_2(grid_cell);
	auto overlap = poly_overlap(*cell_poly, clip_poly);

	return (overlap.size() >= 3);
}

bool EuclidianClip::keep_all(Cell const &grid_cell)
	{ return true; }


}
