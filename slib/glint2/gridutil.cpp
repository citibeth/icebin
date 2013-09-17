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

//#include <snowdrift/maputils.hpp>
#include <glint2/gridutil.hpp>
#include <giss/Proj.hpp>
#include <giss/constant.hpp>
#include <boost/bind.hpp>

namespace glint2 {


Vertex *VertexCache::add_vertex(double x, double y) {
	auto ii = _vertices.find(std::make_pair(x,y));
	if (ii != _vertices.end()) {	// Found in cache
		return ii->second;
	} else {
		Vertex *vertex = grid->add_vertex(Vertex(x,y));
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
