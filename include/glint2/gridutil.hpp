#pragma once

#include <unordered_map>

//#include "geometry.hpp"
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

#include <giss/Proj.hpp>
#include <glint2/Grid.hpp>

// See: http://stackoverflow.com/questions/7222143/unordered-map-hash-function-c
namespace std
{
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

}	// namespace glint2
