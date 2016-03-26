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

#include <unordered_map>

//#include "geometry.hpp"
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

#include <giss/Proj.hpp>
#include <glint2/Grid.hpp>

/**@file std::hash<std::pair<>> is defined in "namespace std" to fix a shortcoming in C++.

@see: http://stackoverflow.com/questions/7222143/unordered-map-hash-function-c
*/
namespace std
{
  /** Used to hash elements of type std::pair<>.  This should have been
  in the standard C++11
  @see: http://stackoverflow.com/questions/7222143/unordered-map-hash-function-c */
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      boost::hash_combine(seed, v.first);
      boost::hash_combine(seed, v.second);
      return seed;
    }
  };
}

namespace glint2 {



/** Tries to eliminate duplicate vertices in the Grid */
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
