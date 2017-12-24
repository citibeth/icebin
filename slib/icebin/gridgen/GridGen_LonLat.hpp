/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ICEBIN_GRID_LONLAT_HPP
#define ICEBIN_GRID_LONLAT_HPP

#include <functional>
#include <memory>
#include <ibmisc/indexing.hpp>

namespace icebin {

/** Represents a Cartesian grid with non-equally-spaced grid cell boundaries. */
extern void Grid make_grid(
    std::string const &name,
    ibmisc::Indexing const &indexing,
    GridSpec_LonLat const &spec,
    std::function<bool(Cell const &)> spherical_clip = &SphericalClip::keep_all);



}   // Namespace
#endif  // Guard
