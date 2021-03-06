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

#pragma once

#include <memory>
#include <functional>
#include <ibmisc/blitz.hpp>
#include <icebin/Grid.hpp>
#include <ibmisc/indexing.hpp>
#include <icebin/gridgen/clippers.hpp>

namespace icebin {

/** Represents a Cartesian grid with non-equally-spaced grid cell boundaries. */
extern Grid make_grid(
    std::string const &name,
    GridSpec_XY const &spec,
    std::function<bool(Cell const &)> const &euclidian_clip = &EuclidianClip::keep_all);



}
