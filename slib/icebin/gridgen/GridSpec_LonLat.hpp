/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
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

#ifndef ICEBIN_GRID_LONLAT_HPP
#define ICEBIN_GRID_LONLAT_HPP

#include <functional>
#include <memory>
#include <ibmisc/indexing.hpp>

namespace icebin {

/** Represents a lat-lon grid on a sphere with non-equally-spaced grid
cell boundaries.  Can optionally represent circular "cap" cells at the
north and south pole, if needed. */
struct GridSpec_LonLat {

    std::string name;
    std::function<bool(double, double, double, double)> spherical_clip;
    ibmisc::Indexing<int,long> indexing;

    /** Longitude of cell boundaries (degrees), sorted low to high.
    <b>NOTE:</b> lon_boundares.last() = 360.0 + lonb.first() */
    std::vector<double> lonb;

    /** Latitude of cell boundaries (degrees), sorted low to high.
    Does NOT include the polar "cap" cells.
    90 = north pole, -90 = south pole. */
    std::vector<double> latb;

    /** True if this grid contains a circular cap on the south pole.
    If not, then many triangular grid cells will meet at the south pole. */
    bool south_pole;

    /** True if this grid contains a circular cap on the north pole.
    If not, then many triangular grid cells will meet at the north pole. */
    bool north_pole;

    /** Number of segments used to represent each side of a grid cell
    with a polygonal approximation.  Increasing this number will
    decrease geometric error at the expense of computation time for
    the overlap matrix. */
    int points_in_side;

    /** Radius of Earth (m), or whatever planet you're looking at.
    Used to compute theoretical exact areas of graticules. That's different
    from the radius (or elliptical shape) used in the projection. */
    double eq_rad;

    // --------------------------------------------

    /** Number of grid cell indices in longitude dimension */
    int nlon() const { return lonb.size() - 1; }

    /** Number of grid cell indices in latitude dimension */
    int nlat() const {
        const int south_pole_offset = (south_pole ? 1 : 0);
        const int north_pole_offset = (north_pole ? 1 : 0);
        return latb.size() - 1 + south_pole_offset + north_pole_offset;
    }

    /** @return [nlat()] latidue of cell centers */
    std::vector<double> latc() const;
    /** @return [nlon()] longitude of cell centers */
    std::vector<double> lonc() const;

    // ------------------------------------------------------

    void make_grid(Grid_LonLat &grid);
};


}   // Namespace
#endif  // Guard
