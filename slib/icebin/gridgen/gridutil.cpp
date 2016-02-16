#include <icebin/gridgen/gridutil.hpp>

namespace icebin {


Vertex *VertexCache::add_vertex(double x, double y) {
	auto ii = _vertices.find(std::make_pair(x,y));
	if (ii != _vertices.end()) {	// Found in cache
		return ii->second;
	} else {
		Vertex *vertex = grid->vertices.add(Vertex(x,y));
		_vertices.insert(std::make_pair(
			std::make_pair(x,y), vertex));
		return vertex;
	}
}

Vertex *VertexCache::add_vertex(Cell &cell, double x, double y)
{
	Vertex *v = add_vertex(x,y);	// Add to the Grid
	cell.add_vertex(v);
	return v;
}


}