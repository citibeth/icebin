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
See Surveyor's Formula: http://www.maa.org/pubs/Calc_articles/ma063.pdf */
extern double area_of_proj_polygon(Cell const &cell, giss::Proj2 const &proj)
{
	double ret = 0;
	auto it0 = cell.begin();
	auto it1(it0); ++it1;

	double x00, y00;
	proj.transform(it0->x, it0->y, x00, y00);

	double x0, y0, x1, y1;
	x0 = x00; y0 = y00;
	for(; it1 != cell.end(); it0=it1, it1 += 1) {
		double x1, y1;
		proj.transform(it1->x, it1->y, x1, y1);

		ret += (x0 * y1) - (x1 * y0);
	}

	x1 = x00;
	y1 = y00;
	ret += (x0 * y1) - (x1 * y0);
	ret *= .5;
	return ret;
}

// ------------------------------------------------------------
Grid::Grid(std::string const &_stype) : stype(_stype) {}

void Grid::clear()
{
	_vertices.clear();
	_cells.clear();
}

Cell *Grid::add_cell(Cell &&cell) {
	// If we never specify our indices, things will "just work"
	if (cell.index == -1) cell.index = _cells.size();


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
void Grid::netcdf_write(NcFile *nc, std::string const &vname) const
{
	// ---------- Write out the vertices
	NcVar *vertices_index_var = nc->get_var((vname + ".vertices.index").c_str());
	NcVar *vertices_xy_var = nc->get_var((vname + ".vertices.xy").c_str());

	std::vector<Vertex *> vertices(_vertices.sorted());
	int i=0;
	for (auto vertex = vertices.begin(); vertex != vertices.end(); ++i, ++vertex) {
		vertices_index_var->set_cur(i);
		vertices_index_var->put(&(*vertex)->index, 1);

		double point[2] = {(*vertex)->x, (*vertex)->y};
		vertices_xy_var->set_cur(i, 0);
		vertices_xy_var->put(point, 1, 2);
	}

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

	std::vector<Cell *> cells(_cells.sorted());
	int ivref = 0;
	i=0;
	for (auto celli = cells.begin(); celli != cells.end(); ++i, ++celli) {
		Cell *cell(*celli);

		// Write general cell contents
		cells_index_var->set_cur(i);
		cells_index_var->put(&cell->index, 1);

// 		cells_i_var->set_cur(i);
// 		cells_i_var->put(&cell->i, 1);
// 
// 		cells_j_var->set_cur(i);
// 		cells_j_var->put(&cell->j, 1);
// 
// 		cells_k_var->set_cur(i);
// 		cells_k_var->put(&cell->k, 1);

		int ijk[3] = {cell->i, cell->j, cell->k};
		cells_ijk_var->set_cur(i, 0);
		cells_ijk_var->put(ijk, 1, 3);

		cells_area_var->set_cur(i);
		cells_area_var->put(&cell->area, 1);

// 		cells_native_area_var->set_cur(i);
// 		double native_area = cell->native_area();
// 		cells_native_area_var->put(&native_area, 1);
// 
// 		cells_proj_area_var->set_cur(i);
// 		double proj_area = cell->proj_area();
// 		cells_proj_area_var->put(&proj_area, 1);

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
}



boost::function<void ()> Grid::netcdf_define(NcFile &nc, std::string const &vname) const
{
printf("netcdf_define(%s)\n", vname.c_str());

	// ------ Attributes
	auto one_dim = giss::get_or_add_dim(nc, "one", 1);
	NcVar *info_var = nc.add_var((vname + ".info").c_str(), ncInt, one_dim);
		info_var->add_att("name", name.c_str());
		info_var->add_att("type", stype.c_str());
		info_var->add_att("coordinates", scoord.c_str());
		if (scoord == "xy") {
			info_var->add_att("projection", sproj.c_str());
//			giss::Proj proj(sproj);
//			giss::Proj llproj(proj.latlong_from_proj());
//			info_var->add_att("llprojection", llproj.get_def().c_str());
		}
		info_var->add_att("cells.num_full", ncells_full);
		info_var->add_att("vertices.num_full", nvertices_full);

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

	return boost::bind(&Grid::netcdf_write, this, &nc, vname);
}

/** @param fname Name of file to load from (eg, an overlap matrix file)
@param vname Eg: "grid1" or "grid2" */
void Grid::read_from_netcdf(
NcFile &nc,
std::string const &vname)
{
	clear();

	printf("Grid::read_from_netcdf(%s)\n", vname.c_str());
	// ---------- Read the Basic Info
	NcVar *info_var = nc.get_var((vname + ".info").c_str());
printf("AA\n");
		name = std::string(info_var->get_att("name")->as_string(0));
		stype = std::string(info_var->get_att("type")->as_string(0));
		scoord = std::string(info_var->get_att("coordinates")->as_string(0));
printf("AA\n");
		if (scoord == "xy")
			sproj = std::string(info_var->get_att("projection")->as_string(0));
		else
			sproj = "";
		ncells_full = info_var->get_att("cells.num_full")->as_int(0);
		nvertices_full = info_var->get_att("vertices.num_full")->as_int(0);
printf("AA\n");

printf("Grid::read_from_netcdf 2");
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
		int index = vertices_index[i];
		double x = vertices_xy[i*2];
		double y = vertices_xy[i*2 + 1];
		add_vertex(Vertex(x, y, index));
	}

//printf("Grid::read_from_netcdf 3\n");
	// ---------- Read the Cells
	std::vector<int> cells_index(giss::read_int_vector(nc, vname + ".cells.index"));

	NcVar *cells_ijk_var = nc.get_var((vname + ".cells.ijk").c_str());
	long ncells = cells_ijk_var->get_dim(0)->size();
	std::vector<int> cells_ijk(ncells*3);
	cells_ijk_var->get(&cells_ijk[0], ncells, 3);

//	std::vector<int> cells_i(giss::read_int_vector(nc, vname + ".cells.i"));
//	std::vector<int> cells_j(giss::read_int_vector(nc, vname + ".cells.j"));
//	std::vector<int> cells_k(giss::read_int_vector(nc, vname + ".cells.k"));
	std::vector<double> cells_area(giss::read_double_vector(nc, vname + ".cells.area"));
//	std::vector<double> cells_native_area(giss::read_double_vector(nc, vname + ".cells.native_area"));
//	std::vector<double> cells_proj_area(giss::read_double_vector(nc, vname + ".cells.proj_area"));

	std::vector<int> vrefs(giss::read_int_vector(nc, vname + ".cells.vertex_refs"));
	std::vector<int> vrefs_start(giss::read_int_vector(nc, vname + ".cells.vertex_refs_start"));

//printf("Grid::read_from_netcdf 4");
	// Assemble into Cells
	for (size_t i=0; i < cells_index.size(); ++i) {
		int index = cells_index[i];

		Cell cell;
		cell.index = cells_index[i];
		cell.i = cells_ijk[i*3 + 0];
		cell.j = cells_ijk[i*3 + 1];
		cell.k = cells_ijk[i*3 + 2];

//		cell.i = cells_i[i];
//		cell.j = cells_j[i];
//		cell.k = cells_k[i];
		cell.area = cells_area[i];
//		cell._native_area = cells_native_area[i];
//		cell._proj_area = cells_proj_area[i];

		// Add the vertices
		cell.reserve(vrefs_start[i+1] - vrefs_start[i]);
		for (int j = vrefs_start[i]; j < vrefs_start[i+1]; ++j)
			cell.add_vertex(get_vertex(vrefs[j]));

		// Add thecell to the grid
		add_cell(std::move(cell));
	}
printf("Grid::read_from_netcdf 5\n");
}

void Grid::to_netcdf(std::string const &fname)
{
	NcFile nc(fname.c_str(), NcFile::Replace);

	// Define stuff in NetCDF file
	printf("Defining netCDF file %s\n", fname.c_str());
	auto gridd = netcdf_define(nc, "grid");

	// Write stuff in NetCDF file
	printf("Writing ot netCDF file\n");
	gridd();

	nc.close();
}

// ---------------------------------------------------

}	// namespace
