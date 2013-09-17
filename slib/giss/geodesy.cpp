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

#include <giss/constant.hpp>

namespace giss {

inline double sqr(double x) { return x*x; }

/** See: http://www.cs.nyu.edu/visual/home/proj/tiger/gisfaq.html
@return Distance (in degrees) */
extern double haversine_distance(
double lon1_deg, double lat1_deg,
double lon2_deg, double lat2_deg)
{
        // Convert inputs to degrees
        double lon1 = lon1_deg * D2R;
        double lat1 = lat1_deg * D2R;
        double lon2 = lon2_deg * D2R;
        double lat2 = lat2_deg * D2R;

        // Apply the Haversine Formula
        double dlon = lon2 - lon1;
        double dlat = lat2 - lat1;
        double a = sqr(sin(dlat/2)) + cos(lat1) * cos(lat2) * sqr(sin(dlon/2));
        double c = 2 * atan2( sqrt(a), sqrt(1-a) );
        return c * R2D;         // Convert to degrees
}

}
