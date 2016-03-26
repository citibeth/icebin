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

#include <icebin/gridgen/clippers.hpp>
#include <ibmisc/geodesy.hpp>

namespace icebin {

// We use lon0...lat1 as a poor man's substute for real CG on a spherical surface

bool SphericalClip::azimuthal(
double center_lon, double center_lat, double clip_distance_deg,
double lon0, double lat0, double lon1, double lat1)
{


//printf("distance=%g, clip_distance=%g\n", ibmisc::haversine_distance(center_lon, center_lat, lon0, lat0), clip_distance_deg);

    // This avoids projection-caused grid cell degeneracy
    // far away from the central area we're interested in
    if (ibmisc::haversine_distance(center_lon, center_lat, lon0, lat0) <= clip_distance_deg) return true;
    if (ibmisc::haversine_distance(center_lon, center_lat, lon1, lat0) <= clip_distance_deg) return true;
    if (ibmisc::haversine_distance(center_lon, center_lat, lon1, lat1) <= clip_distance_deg) return true;
    if (ibmisc::haversine_distance(center_lon, center_lat, lon0, lat1) <= clip_distance_deg) return true;
    return false;
}

/** Tells if one point lies within a longitude range */
inline bool in_lon(double min_lon, double max_lon, double x)
{
    while (x > max_lon) x -= 360.;
    while (x < min_lon) x += 360.;
    bool ret = (x <= max_lon);
//printf("%f %s (%f, %f)\n", x, ret ? "in" : "out", min_lon, max_lon);
    return ret;
}

/** Tells if two longitude ranges intersect. */
inline bool lon_intersect(double min0, double max0, double min1, double max1)
{
    return in_lon(min0, max0, min1) || in_lon(min0, max0, max1)
        || in_lon(min1, max1, min0) || in_lon(min1, max1, max0);
}

bool SphericalClip::lonlat(
double min_lon, double min_lat, double max_lon, double max_lat,
double lon0, double lat0, double lon1, double lat1)
{
//printf("--------------\n");
    if (!lon_intersect(lon0, lon1, min_lon, max_lon)) return false;

    if (lat0 < min_lat && lat1 < min_lat) return false;
    if (lat0 > max_lat && lat1 > max_lat) return false;
    return true;
}

bool SphericalClip::keep_all(
double lon0, double lat0, double lon1, double lat1)
    { return true; }


// =============================================================
/** @param clip_poly Only realize grid cells that intersect with this polygon (on the map) */
bool EuclidianClip::poly(gc::Polygon_2 const &clip_poly,
        Cell const &grid_cell)
{
    auto cell_poly = Cell_to_Polygon_2(grid_cell);
    auto overlap = poly_overlap(*cell_poly, clip_poly);

    return (overlap.size() >= 3);
}

bool EuclidianClip::keep_all(Cell const &grid_cell)
    { return true; }


}
