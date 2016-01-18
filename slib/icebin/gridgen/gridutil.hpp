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

}	// namespace glint2
