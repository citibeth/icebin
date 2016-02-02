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

#include <set>
#include <algorithm>
#include <icebin/Grid.hpp>
#include <ibmisc/netcdf.hpp>
#include <boost/bind.hpp>
//#include <giss/constant.hpp>
#include <icebin/error.hpp>

using namespace ibmisc;
using namespace netCDF;

namespace icebin {

// ========================================================

/** Computes area of the cell's polygon.
For cells in Lat/Lon coordinate, proj is the projection to the plane.
Area of the PROJECTED grid cell is returned.

See Surveyor's`g Formula: http://www.maa.org/pubs/Calc_articles/ma063.pdf */
double Cell::proj_area(
	Proj_LL2XY const *proj)	// OPTIONAL
{
	double ret = 0;
	double x0, y0, x1, y1;

	auto ii(end(-1));		// Last item in the vector
	if (proj) {
		proj->transform(ii->x, ii->y, x0, y0);
	} else {
		x0 = ii->x;
		y0 = ii->y;
	}

	for(ii = begin(); ii != end(); ++ii) {
		double x1, y1;
		if (proj) {
			proj->transform(ii->x, ii->y, x1, y1);
		} else {
			x1 = ii->x;
			y1 = ii->y;
		}

		ret += (x0 * y1) - (x1 * y0);
		x0 = x1;
		y0 = y1;
	}
	ret *= .5;
	return ret;
}

// ========================================================

// ------------------------------------------------------------
Grid::Grid() :
	type(Grid::Type::XY),
	coordinates(Grid::Coordinates::XY),
	_max_realized_cell_index(0),
	_max_realized_vertex_index(0) {}

size_t Grid::ndata() const
{
	if (parameterization == Parameterization::L1)
		return vertices.nfull();
	else
		cells.nfull();
}

void Grid::clear()
{
	vertices.clear();
	cells.clear();
}

// ------------------------------------------------------------
struct CmpVertexXY {
	bool operator()(Vertex const *a, Vertex const *b)
	{
		double diff = a->x - b->x;
		if (diff < 0) return true;
		if (diff > 0) return false;
		return (a->y - b->y) < 0;
	}
};

void sort_renumber_vertices(Grid &grid)
{
	// Construct array of Vertex pointers
	std::vector<Vertex *> vertices;
	for (auto vertex = grid.vertices.begin(); vertex != grid.vertices.end(); ++vertex)
		vertices.push_back(&*vertex);

	// Sort it by x and y!
	std::sort(vertices.begin(), vertices.end(), CmpVertexXY());

	// Renumber vertices
	long i=0;
	for (auto vertex = vertices.begin(); vertex != vertices.end(); ++vertex)
		(*vertex)->index = i++;
}

// ------------------------------------------------------------
void Grid::nc_write(netCDF::NcGroup *nc, std::string const &vname) const
{
	// ---------- Write out the vertices
	{
		std::vector<size_t> startp = {0,0};
		std::vector<size_t> countp = {1,2};
		NcVar vertices_index_v = nc->getVar(vname + ".vertices.index");
		NcVar vertices_xy_v = nc->getVar(vname + ".vertices.xy");

		std::vector<Vertex const *> svertices(vertices.sorted());	// Sort by index
		int i=0;
		for (auto ii = svertices.begin(); ii != svertices.end(); ++i, ++ii) {
			auto vertex(*ii);
			startp[0] = i;
			vertices_index_v.putVar(startp, countp, &vertex->index);

			std::array<double, 2> point = {vertex->x, vertex->y};
			vertices_xy_v.putVar(startp, countp, &point[0]);
		}
	}

	// -------- Write out the cells (and vertex references)
	{
		NcVar cells_index_v = nc->getVar(vname + ".cells.index");
		NcVar cells_ijk_v = nc->getVar(vname + ".cells.ijk");
		NcVar cells_native_area_v = nc->getVar(vname + ".cells.native_area");

		NcVar cells_vertex_refs_v = nc->getVar(vname + ".cells.vertex_refs");
		NcVar cells_vertex_refs_start_v = nc->getVar(vname + ".cells.vertex_refs_start");

		std::vector<Cell const *> scells(cells.sorted());
		std::vector<size_t> startp = {0,0};
		std::vector<size_t> countp = {1,3};
		std::vector<size_t> ivref = {0};
		for (auto ii = scells.begin(); ii != scells.end(); ++ii, ++startp[0]) {
			auto cell(*ii);

			// Write general cell contents
			cells_index_v.putVar(startp, countp, &cell->index);

			std::array<int, 3> ijk = {cell->i, cell->j, cell->k};
			cells_ijk_v.putVar(startp, countp, &ijk[0]);

			cells_native_area_v.putVar(startp, countp, &cell->native_area);

			// Write vertex indices for this cell
			cells_vertex_refs_start_v.putVar(startp, countp, &ivref[0]);
			for (auto vertex = cell->begin(); vertex != cell->end();
				++vertex, ++ivref[0])
			{
				cells_vertex_refs_v.putVar(ivref, countp, &vertex->index);
			}

		}

		// Write out a sentinel for polygon index bounds
		cells_vertex_refs_start_v.putVar(startp, countp, &ivref[0]);
	}
}

/** @param fname Name of file to load from (eg, an overlap matrix file)
@param vname Eg: "gridA" or "gridI" */
void Grid::nc_read(
NcGroup *nc,
std::string const &vname)
{
	clear();

	// ---------- Read the Vertices
	{auto vertices_index(nc_read_blitz
		<long, 1>(nc, vname + ".vertices.index"));
	auto vertices_xy(nc_read_blitz
		<double, 2>(nc, vname + ".vertices.xy"));

		// Assemble into vertices
		for (int i=0; i < vertices_index.extent(0); ++i) {
			long index = vertices_index(i);
			double x = vertices_xy(i,0);
			double y = vertices_xy(i,1);
			vertices.add(Vertex(x, y, index));
		}
	}

	// ---------- Read the Cells
	{
		auto cells_index(nc_read_blitz
			<long, 1>(nc, vname + ".cells.index"));
		auto cells_ijk(nc_read_blitz
			<long, 2>(nc, vname + ".cells.ijk"));
		auto native_area(nc_read_blitz
			<long, 1>(nc, vname + ".cells.native_area"));

		// std::vector<double> cells_area(giss::read_double_vector(nc, vname + ".cells.area"));
		auto vrefs(nc_read_blitz
			<long, 1>(nc, vname + ".cells.vertex_refs"));
		auto vrefs_start(nc_read_blitz
			<long, 1>(nc, vname + ".cells.vertex_refs_start"));


		// Assemble into Cells
		for (int i=0; i < cells_index.size(); ++i) {
			long index = cells_index(i);

			Cell cell;
			cell.index = cells_index(i);
			cell.i = cells_ijk(i,0);
			cell.j = cells_ijk(i,1);
			cell.k = cells_ijk(i,2);
			cell.native_area = native_area(i);

			// Add the vertices
			cell.reserve(vrefs_start(i+1) - vrefs_start(i));
			for (int j = vrefs_start(i); j < vrefs_start(i+1); ++j)
				cell.add_vertex(vertices.at(vrefs(j)));

			// Add the cell to the grid
			cells.add(std::move(cell));
		}
	}
}

std::unique_ptr<Grid> new_grid(Grid::Type type)
{
	switch(type.index()) {
		case Grid::Type::GENERIC :
		case Grid::Type::EXCHANGE :
			return std::unique_ptr<Grid>(new Grid());
		case Grid::Type::XY :
			return std::unique_ptr<Grid>(new Grid_XY());
		case Grid::Type::LONLAT :
			return std::unique_ptr<Grid>(new Grid_LonLat());
		default :
			(*icebin_error)(-1,
				"Unrecognized Grid::Type: %s", type.str());
	}
}


std::unique_ptr<Grid> read_grid(NcIO &ncio, std::string const &vname)
{
	Grid::Type type;
	auto info_v = get_or_add_var(ncio, vname + ".info", netCDF::ncInt64, {});
	get_or_put_att_enum(info_v, ncio.rw, "type", type);

	auto grid(new_grid(type));
	grid->ncio(ncio, vname);

	return grid;
}

void Grid::ncio(NcIO &ncio, std::string const &vname)
{
	// ------ Attributes
	auto info_v = get_or_add_var(ncio, vname + ".info", netCDF::ncInt64, {});

	get_or_put_att(info_v, ncio.rw, "name", name);

	int version = 2;
	get_or_put_att(info_v, ncio.rw, "version", ncInt, &version, 1);
	if (ncio.rw == 'r' && version != 2) {
		(*icebin_error)(-1, "Trying to read version %d, I only know how to read version 2 grids from NetCDF", version);
	}

	get_or_put_att_enum(info_v, ncio.rw, "type", type);
	if (ncio.rw == 'w') info_v.putAtt("type.comment",
		"The overall type of grid, controlling the C++ class used "
		"to represent the grid.  See Grid::Type in slib/icebin/Grid.hpp");

	get_or_put_att_enum(info_v, ncio.rw, "coordinates", coordinates);
	if (ncio.rw == 'w') info_v.putAtt("coordinates.comment",
		"The coordinate system used to represent grid vertices "
		"(See Grid::Coordinates in slib/icebin/Grid.hpp.  May be "
		"either XY or LONLAT (longitude comes before latitude).  Note "
		"that this is different from grid.info.type.  A GENERIC grid, "
		"for example, could be expressed in either XY or LONLAT coordinates.");

	get_or_put_att_enum(info_v, ncio.rw, "parameterization", parameterization);
	if (ncio.rw == 'w') info_v.putAtt("parameterization.comment",
		"Indicates how values are interpolated between grid points "
		"(See Grid::Parameterization in  slib/icebin/Grid.hpp).  Most "
		"finite difference models will use L0, while finite element "
		"models would use L1 or something else.");

	indexing.ncio(ncio, ncInt, vname + ".indexing");

	if (coordinates == Coordinates::XY) {
		get_or_put_att(info_v, ncio.rw, "projection", sproj);
		if (ncio.rw == 'w') info_v.putAtt("projection.comment",
			"If grid.info.coordinates = XY, this indicates the projection "
			"used to convert local XY coordinates to LONLAT coordinates on "
			"the surface of the earth.  See http://trac.osgeo.org/proj/Proj.4 "
			"for format of these strings.");
	}

	get_or_put_att(info_v, ncio.rw, "cells.nfull", ncInt64, &cells._nfull, 1);
	if (ncio.rw == 'w') info_v.putAtt("cells.nfull.comment",
		"The total theoretical number of grid cells (polygons) in this "
		"grid.  Depending on grid.info:parameterization, either cells or "
		"vertices will correspond to the dimensionality of the grid's "
		"vector space.");

	get_or_put_att(info_v, ncio.rw, "vertices.nfull", ncInt64, &vertices._nfull, 1);
	if (ncio.rw == 'w') info_v.putAtt("vertices.nfull.comment",
		"The total theoretical of vertices (of polygons) on this grid.");

	// ------- Dimensions
	if (ncio.rw == 'w') {
		// ----------------- WRITE

		// Count the number of times a vertex (any vertex) is referenced.
		int nvref = 0;
		for (auto cell = cells.begin(); cell != cells.end(); ++cell) nvref += cell->size();


		NcDim vertices_nrealized_d = get_or_add_dim(ncio, vname + ".vertices.nrealized", vertices.nrealized());
		info_v.putAtt((vname + ".vertices.nrealized.comment"),
			"The number of 'realized' cells in this grid.  Only the "
			"outlines of realized cells are computed and stored.  not "
			"all cells need to be realized.  For example, a grid file "
			"representing a GCM grid, in preparation for use with ice "
			"models, would only need to realize GCM grid cells that are "
			"close to the relevant ice sheets.  In this case, all grid "
			"cells are realized.");

		NcDim cells_nrealized_d = get_or_add_dim(ncio,
			vname + ".cells.nrealized", cells.nrealized());
		NcDim cells_nrealized_plus_1_d = get_or_add_dim(ncio,
			vname + ".cells.nrealized_plus1", cells.nrealized() + 1);

		NcDim nvrefs_d = get_or_add_dim(ncio, vname + ".cells.nvertex_refs", nvref);
		NcDim two_d = get_or_add_dim(ncio, "two", 2);
		NcDim three_d = get_or_add_dim(ncio, "three", 3);

		// --------- Variables
		get_or_add_var(ncio, vname + ".vertices.index", ncInt, {vertices_nrealized_d})
			.putAtt("comment",
				"For grids that index on cells (eg, L0): a dense, zero-based "
				"1D index used to identify each realized cell.  This will be "
				"used for vectors representing fields on the grid.");

		get_or_add_var(ncio, vname + ".vertices.xy", ncDouble, {vertices_nrealized_d, two_d});

		get_or_add_var(ncio, vname + ".cells.index", ncInt, {cells_nrealized_d})
			.putAtt("comment",
				"For grids that index on vertices (eg, L1): a dense, zero-based "
				"1D index used to identify each realized vertex.  This will be "
				"used for vectors representing fields on the grid.");

		get_or_add_var(ncio, vname + ".cells.ijk", ncInt, {cells_nrealized_d, three_d})
			.putAtt("comment",
				"OPTIONAL: Up to 3 dimensions can be used to assign a 'real-world' "
				"index to each grid cell.  If grid.info:type = EXCHANGE, then i and "
				"j correspond to grid.vertices.index of the two overlapping source cells.");

		get_or_add_var(ncio, vname + ".cells.native_area", ncDouble, {cells_nrealized_d})
			.putAtt("comment",
				"Area of each cell in its native (non-projected) coordinate system.  "
				"We can compute the projected area on the fly.");


		// nc.add_var((vname + ".cells.area").c_str(), ncDouble, ncells_dim);

		get_or_add_var(ncio, vname + ".cells.vertex_refs", ncInt, {nvrefs_d});
		get_or_add_var(ncio, vname + ".cells.vertex_refs_start", ncInt, {cells_nrealized_plus_1_d});

		ncio += std::bind(&Grid::nc_write, this, ncio.nc, vname);
	} else {
		// ----------------- READ
		nc_read(ncio.nc, vname);
	}

}

// ============================================================

/** Creates a new grid with just the cells we want to keep in it.
This is used to remove cells that don't fit our MPI domain. */
void Grid::filter_cells(std::function<bool(long)> const &keep_fn)
{
	std::set<int> good_vertices;	// Remove vertices that do NOT end up in this set.

printf("BEGIN filter_cells(%s) %p\n", name.c_str(), this);

	// Set counts so they won't change
	cells._nfull = cells.nfull();
	vertices._nfull = vertices.nfull();

	// Remove cells that don't fit our filter
	_max_realized_cell_index = -1;
	for (auto cell = cells.begin(); cell != cells.end(); ) { //++cell) {
		if (keep_fn(cell->index)) {
			_max_realized_cell_index = std::max(_max_realized_cell_index, cell->index);

			// Make sure we don't delete this cell's vertices
			for (auto vertex = cell->begin(); vertex != cell->end(); ++vertex)
				good_vertices.insert(vertex->index);
			++cell;
		} else {
			// Remove the cell, maybe remove its vertices later
			// Careful with iterators: invalidated after erase()
			cell = cells.erase(cell);	// Increments too
		}
	}

	// Remove vertices that don't fit our filter
	_max_realized_vertex_index = -1;
	for (auto vertex = vertices.begin(); vertex != vertices.end(); ) {
		if (good_vertices.find(vertex->index) != good_vertices.end()) {
			_max_realized_vertex_index = std::max(_max_realized_vertex_index, vertex->index);
			++vertex;
		} else {
			// Careful with iterators: invalidated after erase()
			vertex = vertices.erase(vertex);	// Increments too
		}
	}

	printf("END filter_cells(%s) %p\n", name.c_str(), this);
}
// ---------------------------------------------------------
// ---------------------------------------------------------
void Grid_XY::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{

	auto xb_d = get_or_add_dim(ncio,
		vname + ".x_boundaries.length", this->xb.size());
	ncio_vector(ncio, this->xb, true,
		vname + ".x_boundaries", ncDouble, {xb_d});

	auto yb_d = get_or_add_dim(ncio,
		vname + ".y_boundaries.length", this->yb.size());
	ncio_vector(ncio, this->yb, true,
		vname + ".y_boundaries", ncDouble, {yb_d});

	NcVar info_v = get_or_add_var(ncio, vname + ".info", ncInt, {});
	if (ncio.rw == 'w') {
		int n;
		n = nx();
		get_or_put_att(info_v, ncio.rw, "nx", ncInt, &n, 1);
		n = ny();
		get_or_put_att(info_v, ncio.rw, "ny", ncInt, &n, 1);
	}

	Grid::ncio(ncio, vname);
}
// ---------------------------------------------------------
void Grid_LonLat::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{

	auto lonb_d = get_or_add_dim(ncio,
		vname + ".lon_boundaries.length", this->lonb.size());
	ncio_vector(ncio, this->lonb, true,
		vname + ".lon_boundaries", ncDouble, {lonb_d});

	auto latb_d = get_or_add_dim(ncio,
		vname + ".lat_boundaries.length", this->latb.size());
	ncio_vector(ncio, this->latb, true,
		vname + ".lat_boundaries", ncDouble, {latb_d});

	NcVar info_v = get_or_add_var(ncio, vname + ".info", ncInt, {});
	get_or_put_att(info_v, ncio.rw, "north_pole_cap", north_pole);
	get_or_put_att(info_v, ncio.rw, "south_pole_cap", south_pole);
	get_or_put_att(info_v, ncio.rw, "points_in_side", ncInt, &points_in_side, 1);
	if (ncio.rw == 'w') {
		int n;
		n = nlon();
		get_or_put_att(info_v, ncio.rw, "nlon", ncInt, &n, 1);
		n = nlat();
		get_or_put_att(info_v, ncio.rw, "nlat", ncInt, &n, 1);
	}

	Grid::ncio(ncio, vname);
}
// ---------------------------------------------------------

int Grid_LonLat::nlat() const {
	int south_pole_offset, north_pole_offset;

printf("poles = %d %d\n", south_pole, north_pole);

	// Get around GCC bug when converting uninitialized bools
	south_pole_offset = (south_pole ? 2 : 0);
	north_pole_offset = (north_pole ? 2 : 0);

	south_pole_offset >>= 1;
	north_pole_offset >>= 1;

printf("Offsets = %d %d\n", south_pole_offset, north_pole_offset);
	return latb.size() - 1 + south_pole_offset + north_pole_offset;
}





}	// namespace

std::ostream &operator<<(std::ostream &os, icebin::Vertex const &vertex)
{
	os << vertex.index << ":(" << vertex.x << ", " << vertex.y << ")"; 
	return os;
}

std::ostream &operator<<(std::ostream &os, icebin::Cell const &cell)
{
	os << "Cell(ix=" << cell.index << ": [";
	for (auto ii(cell.begin()); ii != cell.end(); ++ii) {
		icebin::Vertex const &vertex(*ii);
		os << vertex;
		os << ", ";
	}
	os << "])";
	return os;
}

// ----------------------------------------------------
