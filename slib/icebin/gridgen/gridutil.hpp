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

#pragma once

#include <unordered_map>
#include <functional>
#include <ibmisc/hash.hpp>
#include <icebin/Grid.hpp>

namespace icebin {

/** Eliminates duplicate vertices in the Grid.
DOES NOT account for floating point rounding error */
class VertexCache {
    std::unordered_map<std::pair<double,double>, Vertex *> _vertices;
public:
    Grid *grid;
    VertexCache(Grid *_grid) : grid(_grid) {}

    Vertex *add_vertex(double x, double y);
    Vertex *add_vertex(Cell &cell, double x, double y);
};




// Normalises a value of longitude to the range starting at min degrees.
// @return The normalised value of longitude.
inline double loncorrect(double lon, double min)
{
    double max = min + 360.0;

    while (lon >= max) lon -= 360.0;
    while (lon < min) lon += 360.0;

    return lon;
}

}   // namespace glint2
