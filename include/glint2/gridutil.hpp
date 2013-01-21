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


/** A helper class */
class LineProjector {

public :
	giss::Proj const proj;
	giss::Proj const llproj;
private:
	/** Function to add a vertex to the grid (but not a cell) */
	boost::function<Vertex *(double, double)> _add_vertex;

public :
	LineProjector(giss::Proj const &_proj, giss::Proj const &_llproj,
		boost::function<Vertex *(double, double)> const &add_vertex);

	/** Convenience constructor. */
	LineProjector(giss::Proj const &_proj, VertexCache *vcache);

	/** Project a latitude line segment
	n >=1, number of points to add to polygoon.  Don't add (lat,lon1). */
	void proj_latitude(Cell &cell, int n,
		double lon0, double lon1, double lat);


	/** Project a latitude line segment
	n >=1, number of points to add to polygoon.  Don't add (lat,lon1). */
	void proj_meridian(Cell &cell, int n,
		double lon, double lat0, double lat1);
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
