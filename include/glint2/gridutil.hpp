#pragma once

#include "geometry.hpp"
#include "Proj.hpp"

namespace glint2 {



/** Tries to eliminate duplicate vertices in the Grid */
class VertexCache {
	std::unordered_map<std::pair<double,double>, Vertex *> _vertices;
public:
	Grid const *grid;
	VertexCache(Grid *_grid) : grid(_grid) {}

	Vertex *add_vertex(double x, double y) {
		auto ii = _vertices.find(std::make_pair(x,y));
		if (ii != _vertices.end()) {	// Found in cache
			return ii->second;
		} else {
			return grid->add_vertex(Vertex(x,y));
		}
	}

	Vertex *add_vertex(Cell &cell, double x, double y)
	{
		Vertex *v = add_vertex(x,y);	// Add to the Grid
		cell.add_vertex(v);
		return v;
	}

};



/** A helper class */
class LineProjector {

public :
	Proj const &proj;
	Proj const &llproj;
private:
	boost::function<double, double> _add_vertex;

public :
	LineProjector(Proj const &_proj, Proj const &_llproj,
		boost::function<double, double> _add_vertex) :
		proj(_proj), llproj(_llproj), _add_vertex(add_vertex)
	{}

	/** Convenience constructor. */
	LineProjector(Proj const &_proj, VertexCache *vcache) :
		proj(_proj),
		llproj(proj.latlong_from_proj()),
		_add_vertex(boost::bind(&VertexCache::add_vertex, vcache, _1, _2))
	{}

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
