#include <set>
#include <algorithm>
#include <glint2/Grid.hpp>
#include <giss/ncutil.hpp>
#include <boost/bind.hpp>
#include <giss/constant.hpp>

namespace glint2 {



/** See Surveyor's Formula: http://www.maa.org/pubs/Calc_articles/ma063.pdf */
extern double area_of_polygon(Cell const &cell)
{
	double ret = 0;
	auto it0 = cell.begin();
	auto it1(it0); ++it1;
	for(; it1 != cell.end(); it0=it1, it1 += 1) {
		ret += (it0->x * it1->y) - (it1->x * it0->y);
	}
	it1 = cell.begin();
	ret += (it0->x * it1->y) - (it1->x * it0->y);
	ret *= .5;
	return ret;
}

/** Computes area of the cell's polygon after it's projected
(for cells in lon/lat coordinates)
See Surveyor's`g Formula: http://www.maa.org/pubs/Calc_articles/ma063.pdf */
extern double area_of_proj_polygon(Cell const &cell, giss::Proj2 const &proj)
{
	double ret = 0;
	double x00, y00;
	auto it0 = cell.begin();
	proj.transform(it0->x, it0->y, x00, y00);

//printf("(%f, %f) --> (%f, %f)\n", it0->x, it0->y, x00, y00);

	double x0, y0, x1, y1;
	x0 = x00; y0 = y00;
	for(++it0; it0 != cell.end(); ++it0) {
		double x1, y1;
		proj.transform(it0->x, it0->y, x1, y1);

		ret += (x0 * y1) - (x1 * y0);
		x0 = x1;
		y0 = y1;
	}

	x1 = x00;
	y1 = y00;
	ret += (x0 * y1) - (x1 * y0);
	ret *= .5;
	return ret;
}

/** Compute the center of a polygon in a plane.
NOTE: Does NOT work (correctly) for lon/lat coordinates
http://stackoverflow.com/questions/5271583/center-of-gravity-of-a-polygon */
extern void center_of_polygon(Cell const &cell, double &x, double &y)
{
	double A = area_of_polygon(cell);
	double Cx = 0;
	double Cy = 0;
	auto v0(cell.begin());
	auto v1(v0);
	for (++ v1; v1 != cell.end(); v0=v1,++v1) {
		double cross = (v0->x * v1->y - v1->x * v0->y);
		Cx += (v0->x + v1->x) * cross;
		Cy += (v0->y + v1->y) * cross;
	}

	double fact = 1.0 / (6.0 * A);
	x= Cx * fact;
	y = Cy * fact;
}

// ------------------------------------------------------------
Grid::Grid(Type _type) :
	type(_type),
	coordinates(Grid::Coordinates::LONLAT),
	parameterization(Grid::Parameterization::L0),
	_max_realized_cell_index(0),
	_max_realized_vertex_index(0) {}

long Grid::ndata() const
{
	if (parameterization == Parameterization::L1) {
		return nvertices_full();
	}
	return ncells_full();
}

void Grid::clear()
{
	_vertices.clear();
	_cells.clear();
}

Cell *Grid::add_cell(Cell &&cell) {
	// If we never specify our indices, things will "just work"
	if (cell.index == -1) cell.index = _cells.size();
	_max_realized_cell_index = std::max(_max_realized_cell_index, cell.index);

	auto ret = _cells.insert(cell.index, std::move(cell));
	Cell *valp = ret.first;
	bool inserted = ret.second;

	if (!inserted) {		// Key already existed
		fprintf(stderr, "Error adding repeat cell index=%d.  "
			"Cells must have unique indices.", cell.index);
		throw std::exception();
	}
	return valp;
}

Vertex *Grid::add_vertex(Vertex &&vertex) {
	// If we never specify our indices, things will "just work"
	if (vertex.index == -1) vertex.index = _vertices.size();
	_max_realized_vertex_index = std::max(_max_realized_vertex_index, vertex.index);

//printf("add_vertex(%d, %f, %f)\n", vertex.index, vertex.x, vertex.y);

	auto ret = _vertices.insert(vertex.index, std::move(vertex));
	Vertex *valp = ret.first;
	bool inserted = ret.second;

	if (!inserted) {		// Key already existed
		fprintf(stderr, "Error adding repeat vertex index=%d.  "
			"Vertices must have unique indices.", vertex.index);
		throw std::exception();
	}
	return valp;
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

void Grid::sort_renumber_vertices()
{
	// Construct array of Vertex pointers
	std::vector<Vertex *> vertices;
	for (auto vertex = vertices_begin(); vertex != vertices_end(); ++vertex)
		vertices.push_back(&*vertex);

	// Sort it by x and y!
	std::sort(vertices.begin(), vertices.end(), CmpVertexXY());

	// Renumber vertices
	long i=0;
	for (auto vertex = vertices.begin(); vertex != vertices.end(); ++vertex)
		(*vertex)->index = i++;
}

// ------------------------------------------------------------
void Grid::netcdf_write(NcFile *nc, std::string const &vname) const
{
printf("Grid::netcdf_write(%s) 1\n", vname.c_str());
	// ---------- Write out the vertices
	NcVar *vertices_index_var = nc->get_var((vname + ".vertices.index").c_str());
	NcVar *vertices_xy_var = nc->get_var((vname + ".vertices.xy").c_str());

	std::vector<Vertex *> vertices(_vertices.sorted());	// Sort by index
	int i=0;
	for (auto vertex = vertices.begin(); vertex != vertices.end(); ++i, ++vertex) {
		vertices_index_var->set_cur(i);
		vertices_index_var->put(&(*vertex)->index, 1);

		double point[2] = {(*vertex)->x, (*vertex)->y};
		vertices_xy_var->set_cur(i, 0);
		vertices_xy_var->put(point, 1, 2);
	}

printf("Grid::netcdf_write() 2\n");
	// -------- Write out the cells (and vertex references)
	NcVar *cells_index_var = nc->get_var((vname + ".cells.index").c_str());
//	NcVar *cells_i_var = nc->get_var((vname + ".cells.i").c_str());
//	NcVar *cells_j_var = nc->get_var((vname + ".cells.j").c_str());
//	NcVar *cells_k_var = nc->get_var((vname + ".cells.k").c_str());
	NcVar *cells_ijk_var = nc->get_var((vname + ".cells.ijk").c_str());
	NcVar *cells_area_var = nc->get_var((vname + ".cells.area").c_str());
//	NcVar *cells_native_area_var = nc->get_var((vname + ".cells.native_area").c_str());
//	NcVar *cells_proj_area_var = nc->get_var((vname + ".cells.proj_area").c_str());

	NcVar *cells_vertex_refs_var = nc->get_var((vname + ".cells.vertex_refs").c_str());
	NcVar *cells_vertex_refs_start_var = nc->get_var((vname + ".cells.vertex_refs_start").c_str());

printf("Grid::netcdf_write() 3\n");
	std::vector<Cell *> cells(_cells.sorted());
	int ivref = 0;
	i=0;
	for (auto celli = cells.begin(); celli != cells.end(); ++i, ++celli) {
		Cell *cell(*celli);

		// Write general cell contents
		cells_index_var->set_cur(i);
		cells_index_var->put(&cell->index, 1);

		int ijk[3] = {cell->i, cell->j, cell->k};
		cells_ijk_var->set_cur(i, 0);
		cells_ijk_var->put(ijk, 1, 3);

		cells_area_var->set_cur(i);
		cells_area_var->put(&cell->area, 1);

		// Write vertex indices for this cell
		cells_vertex_refs_start_var->set_cur(i);
		cells_vertex_refs_start_var->put(&ivref, 1);
		for (auto vertex = cell->begin(); vertex != cell->end(); ++vertex) {
			cells_vertex_refs_var->set_cur(ivref);
			cells_vertex_refs_var->put(&vertex->index, 1);
			++ivref;
		}
	}

	// Write out a sentinel for polygon index bounds
	cells_vertex_refs_start_var->set_cur(i);
	cells_vertex_refs_start_var->put(&ivref, 1);
printf("Grid::netcdf_write() 5\n");
}



boost::function<void ()> Grid::netcdf_define(NcFile &nc, std::string const &vname) const
{
printf("netcdf_define(%s) 1\n", vname.c_str());

	// ------ Attributes
	auto one_dim = giss::get_or_add_dim(nc, "one", 1);
	NcVar *info_var = nc.add_var((vname + ".info").c_str(), ncInt, one_dim);
		info_var->add_att("version", 1);		// File format versioning
		info_var->add_att("name", name.c_str());
		info_var->add_att("type", type.str());
		info_var->add_att("coordinates", coordinates.str());
		info_var->add_att("parameterization", parameterization.str());
		if (coordinates == Coordinates::XY) {
			info_var->add_att("projection", sproj.c_str());
//			giss::Proj proj(sproj);
//			giss::Proj llproj(proj.latlong_from_proj());
//			info_var->add_att("llprojection", llproj.get_def().c_str());
		}

		char buf[32];
		sprintf(buf, "%ld", ncells_full());
		info_var->add_att("cells.num_full", buf);
		info_var->add_att("vertices.num_full", (int)nvertices_full());
printf("CC\n");

printf("netcdf_define(%s) 2\n", vname.c_str());
	// ------- Dimensions
	// Count the number of times a vertex (any vertex) is referenced.
	int nvref = 0;
	for (auto cell = cells_begin(); cell != cells_end(); ++cell) {
//printf("Found cell %d (size=%d)\n", cell->index, cell->size());
		nvref += cell->size();
	}

	NcDim *nvertices_dim = nc.add_dim(
		(vname + ".vertices.num_realized").c_str(), nvertices_realized());
	NcDim *ncells_dim = nc.add_dim(
		(vname + ".cells.num_realized").c_str(), ncells_realized());
	NcDim *ncells_plus_1_dim = nc.add_dim(
		(vname + ".cells.num_realized_plus1").c_str(), ncells_realized() + 1);
	NcDim *nvrefs_dim = nc.add_dim(
		(vname + ".cells.num_vertex_refs").c_str(), nvref);
	NcDim *two_dim = giss::get_or_add_dim(nc, "two", 2);
	NcDim *three_dim = giss::get_or_add_dim(nc, "three", 3);

printf("netcdf_define(%s) 3\n", vname.c_str());
	// --------- Variables
	nc.add_var((vname + ".vertices.index").c_str(), ncInt, nvertices_dim);
	nc.add_var((vname + ".vertices.xy").c_str(), ncDouble, nvertices_dim, two_dim);

	nc.add_var((vname + ".cells.index").c_str(), ncInt, ncells_dim);
	nc.add_var((vname + ".cells.ijk").c_str(), ncInt, ncells_dim, three_dim);
//	nc.add_var((vname + ".cells.i").c_str(), ncInt, ncells_dim);
//	nc.add_var((vname + ".cells.j").c_str(), ncInt, ncells_dim);
//	nc.add_var((vname + ".cells.k").c_str(), ncInt, ncells_dim);
	nc.add_var((vname + ".cells.area").c_str(), ncDouble, ncells_dim);
//	nc.add_var((vname + ".cells.native_area").c_str(), ncDouble, ncells_dim);
//	nc.add_var((vname + ".cells.proj_area").c_str(), ncDouble, ncells_dim);

	nc.add_var((vname + ".cells.vertex_refs").c_str(), ncInt, nvrefs_dim);
	nc.add_var((vname + ".cells.vertex_refs_start").c_str(), ncInt, ncells_plus_1_dim);

printf("netcdf_define(%s) 4\n", vname.c_str());
	return boost::bind(&Grid::netcdf_write, this, &nc, vname);
}

/** @param fname Name of file to load from (eg, an overlap matrix file)
@param vname Eg: "grid1" or "grid2" */
void Grid::read_from_netcdf(
NcFile &nc,
std::string const &vname)
{
	clear();

	printf("Grid::read_from_netcdf(%s) 1\n", vname.c_str());
	// ---------- Read the Basic Info
	NcVar *info_var = nc.get_var((vname + ".info").c_str());
		name = std::string(giss::get_att(info_var, "name")->as_string(0));

		type = giss::parse_enum<decltype(type)>(
			giss::get_att(info_var, "type")->as_string(0));
		coordinates = giss::parse_enum<decltype(coordinates)>(
			giss::get_att(info_var, "coordinates")->as_string(0));
		parameterization = giss::parse_enum<decltype(parameterization)>(
			giss::get_att(info_var, "parameterization")->as_string(0));

		if (coordinates == Coordinates::XY)
			sproj = std::string(giss::get_att(info_var, "projection")->as_string(0));
		else
			sproj = "";

		char *sncells_full = giss::get_att(info_var, "cells.num_full")->as_string(0);
		sscanf(sncells_full, "%ld", &_ncells_full);

		_nvertices_full = giss::get_att(info_var, "vertices.num_full")->as_int(0);

	printf("Grid::read_from_netcdf(%s) 2\n", vname.c_str());

	// ---------- Read the Vertices
	// Basic Info
	std::vector<int> vertices_index(
		giss::read_int_vector(nc, vname + ".vertices.index"));

	// Read points 2-d array as single vector (double)
	NcVar *vertices_xy_var = nc.get_var((vname + ".vertices.xy").c_str());
	long npoints = vertices_xy_var->get_dim(0)->size();
	std::vector<double> vertices_xy(npoints*2);
	vertices_xy_var->get(&vertices_xy[0], npoints, 2);

	// Assemble into vertices
	for (size_t i=0; i < vertices_index.size(); ++i) {
		long index = vertices_index[i];
		double x = vertices_xy[i*2];
		double y = vertices_xy[i*2 + 1];
		add_vertex(Vertex(x, y, index));
	}
	printf("Grid::read_from_netcdf(%s) 3\n", vname.c_str());

	// ---------- Read the Cells
	std::vector<int> cells_index(giss::read_int_vector(nc, vname + ".cells.index"));

	NcVar *cells_ijk_var = nc.get_var((vname + ".cells.ijk").c_str());
	long ncells = cells_ijk_var->get_dim(0)->size();
	std::vector<int> cells_ijk(ncells*3);
	cells_ijk_var->get(&cells_ijk[0], ncells, 3);

	std::vector<double> cells_area(giss::read_double_vector(nc, vname + ".cells.area"));

	std::vector<int> vrefs(giss::read_int_vector(nc, vname + ".cells.vertex_refs"));
	std::vector<int> vrefs_start(giss::read_int_vector(nc, vname + ".cells.vertex_refs_start"));

	printf("Grid::read_from_netcdf(%s) 4\n", vname.c_str());

	// Assemble into Cells
	for (size_t i=0; i < cells_index.size(); ++i) {
		long index = cells_index[i];

		Cell cell;
		cell.index = cells_index[i];
		cell.i = cells_ijk[i*3 + 0];
		cell.j = cells_ijk[i*3 + 1];
		cell.k = cells_ijk[i*3 + 2];

		cell.area = cells_area[i];

		// Add the vertices
		cell.reserve(vrefs_start[i+1] - vrefs_start[i]);
		for (int j = vrefs_start[i]; j < vrefs_start[i+1]; ++j)
			cell.add_vertex(get_vertex(vrefs[j]));

		// Add thecell to the grid
		add_cell(std::move(cell));
	}
	printf("Grid::read_from_netcdf(%s) 5 (done)\n", vname.c_str());
}

void Grid::to_netcdf(std::string const &fname)
{
	NcFile nc(fname.c_str(), NcFile::Replace);

	// Define stuff in NetCDF file
	printf("Defining netCDF file %s\n", fname.c_str());
	auto gridd = netcdf_define(nc, "grid");

	// Write stuff in NetCDF file
	printf("Writing to netCDF file: %s\n", fname.c_str());
	gridd();

	nc.close();
}

// ---------------------------------------------------
static double const nan = std::numeric_limits<double>::quiet_NaN();

std::vector<double> Grid::get_native_areas() const
{
	// Get the cell areas
	std::vector<double> area(this->ncells_full(), nan);
	for (auto cell = this->cells_begin(); cell != this->cells_end(); ++cell) {
		area.at(cell->index) = cell->area;
	}

	return area;
}

void Grid::get_ll_to_xy(giss::Proj2 &proj, std::string const &sproj) const
{
printf("get_ll_to_xy(sproj=%s)\n", sproj.c_str());
	// Set up the projection
	if (coordinates == Coordinates::LONLAT) {
		proj.init(sproj, giss::Proj2::Direction::LL2XY);
	} else {
		fprintf(stderr, "get_ll_to_xy() only makes sense for grids in Lon/Lat Coordinates!");
		throw std::exception();
	}
}


std::vector<double> Grid::get_proj_areas(std::string const &sproj) const
{
	giss::Proj2 proj;
	get_ll_to_xy(proj, sproj);

	// Get the projected cell areas
	std::vector<double> area(this->ncells_full(), nan);
	for (auto cell = this->cells_begin(); cell != this->cells_end(); ++cell) {
		area.at(cell->index) = area_of_proj_polygon(*cell, proj);
	}

	return area;
}
// ---------------------------------------------------
/** Remove cells and vertices not relevant to us --- for example, not in our MPI domain. */
void Grid::filter_cells(boost::function<bool (int)> const &include_cell)
{
	std::set<int> good_vertices;	// Remove vertices that do NOT end up in this set.

printf("BEGIN filter_cells(%s) %p\n", name.c_str(), this);

	// Set counts so they won't change
	_ncells_full = ncells_full();
	_nvertices_full = nvertices_full();

	// Remove cells that don't fit our filter
	_max_realized_cell_index = -1;
	for (auto cell = cells_begin(); cell != cells_end(); ) { //++cell) {
		bool keep = include_cell(cell->index);
		if (keep) {
			_max_realized_cell_index = std::max(_max_realized_cell_index, cell->index);

			// Make sure we don't delete this cell's vertices
			for (auto vertex = cell->begin(); vertex != cell->end(); ++vertex)
				good_vertices.insert(vertex->index);
			++cell;
		} else {
			// Remove the cell, maybe remove its vertices later
			// Careful with iterators: invalidated after erase()
			cell = _cells.erase(cell);	// Increments too
		}
	}

	// Remove vertices that don't fit our filter
	_max_realized_vertex_index = -1;
	for (auto vertex = vertices_begin(); vertex != vertices_end(); ) {
		if (good_vertices.find(vertex->index) != good_vertices.end()) {
			_max_realized_vertex_index = std::max(_max_realized_vertex_index, vertex->index);
			++vertex;
		} else {
			// Careful with iterators: invalidated after erase()
			vertex = _vertices.erase(vertex);	// Increments too
		}
	}

	printf("END filter_cells(%s) %p\n", name.c_str(), this);
}

}	// namespace
