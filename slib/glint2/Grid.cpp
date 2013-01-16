namespace glint2 {


static void Grid_netcdf_write(Grid const &grid, NcFile *nc, std::string const &vname) const
{
	// ---------- Write out the vertices
	NcVar *vertices_index_var = nc->get_var((vname + ".vertices.index").c_str());
	NcVar *vertices_xy_var = nc->get_var((vname + ".vertices.xy").c_str());

	std::vector<Vertex *> vertices(_vertices.sorted());
	int i=0;
	for (int i=0, auto ii = _vertices.begin(); ii != _vertices.end(); ++i, ++ii) {
		Vertex const *v = *ii;

		vertices_index_var->set_cur(i);
		vertices_index_var->put(&v->index, 1);

		double point[2] = {v->x, v->y};
		vertices_xy_var->set_cur(i, 0);
		vertices_xy_var->put(point, 1, 2);
	}

	// -------- Write out the cells (and vertex references)
	NcVar *cells_index_var = nc->get_var((vname + ".cells.index").c_str());
	NcVar *cells_i_var = nc->get_var((vname + ".cells.i").c_str());
	NcVar *cells_j_var = nc->get_var((vname + ".cells.j").c_str());
	NcVar *cells_k_var = nc->get_var((vname + ".cells.k").c_str());
	NcVar *cells_native_area_var = nc->get_var((vname + ".cells.native_area").c_str());
	NcVar *cells_proj_area_var = nc->get_var((vname + ".cells.proj_area").c_str());

	NcVar *cells_vertex_refs_var = nc->get_var((vname + ".cells.vertex_refs").c_str());
	NcVar *cells_vertex_refs_start_var = nc->get_var((vname + ".cells.vertex_refs_start").c_str());

	int ivref = 0;
	i=0;
	for (auto ii = _cells.begin(); ii != _cells.end(); ++i, ++ii) {
		Cell const *cell = *ii;

		// Write general cell contents
		cells_index_var->set_cur(i);
		cells_index_var->put(&cell->index, 1);

		cells_i_var->set_cur(i);
		cells_i_var->put(&cell->i, 1);

		cells_j_var->set_cur(i);
		cells_j_var->put(&cell->j, 1);

		cells_k_var->set_cur(i);
		cells_k_var->put(&cell->k, 1);

		cells_native_area_var->set_cur(i);
		cells_native_area_var->put(&cell->native_area, 1);

		cells_proj_area_var->set_cur(i);
		cells_proj_area_var->put(&cell->proj_area, 1);

		// Write vertex indices for this cell
		cells_vertex_refs_start_var->set_cur(i);
		cells_vertex_refs_start_var->put(&ivref, 1);
		for (auto vertex = cell->begin(); vertex != cell->end(); ++vertex) {
			cells_vertex_refs_var->set_cur(ivref);
			cells_vertex_refs_var->put(&vertex->index, 1);
		}
	}

	// Write out a sentinel for polygon index bounds
	cells_vertex_refs_start_var->set_cur(i);
	cells_vertex_refs_start_var->put(&ivref, 1);
}



boost::function<void ()> Grid::netcdf_define(NcFile &nc, std::string const &vname) const
{
	// ------ Attributes
	auto one_dim = get_or_add_dim(nc, "one", 1);
	NcVar *info_var = nc.add_var((vname + ".info").c_str(), ncInt, one_dim);
		info_var->add_att("name", name.c_str());
		info_var->add_att("type", stype.c_str());
		info_var->add_att("cells.num_full", ncells());
		info_var->add_att("vertices.num_full", nvertices());

	// ------- Dimensions
	// Count the number of times a vertex (any vertex) is referenced.
	int nvref = 0;
	for (auto ii = cells_begin(); ii != cells_end(); ++ii)
		nvref += ii->size();

	NcDim *nvertices_dim = nc.add_dim((vname + ".vertices.num_realized").c_str(), nvertices_realized());
	NcDim *ncells_dim = nc.add_dim((vname + ".cells.num_realized").c_str(), ncells_realized());
	NcDim *ncells_plus_1_dim = nc.add_dim((vname + ".cells.num_realized_plus1").c_str(), ncells() + 1);
	NcDim *nvrefs_dim = nc.add_dim((vname + ".cells.num_vertex_refs").c_str(), nvref);
	NcDim *two_dim = get_or_add_dim(nc, "two", 2);

	// --------- Variables
	nc.add_var((vname + ".vertices.index").c_str(), ncInt, nvertices_dim);
	nc.add_var((vname + ".vertices.xy").c_str(), ncDouble, nvertices_dim, two_dim);

	nc.add_var((vname + ".cells.index").c_str(), ncInt, ncells_dim);
	nc.add_var((vname + ".cells.i").c_str(), ncInt, ncells_dim);
	nc.add_var((vname + ".cells.j").c_str(), ncInt, ncells_dim);
	nc.add_var((vname + ".cells.k").c_str(), ncInt, ncells_dim);
	nc.add_var((vname + ".cells.native_area").c_str(), ncDouble, ncells_dim);
	nc.add_var((vname + ".cells.proj_area").c_str(), ncDouble, ncells_dim);

	nc.add_var((vname + ".cells.vertex_refs").c_str(), ncInt, nvrefs_dim);
	nc.add_var((vname + ".cells.vertex_refs_start").c_str(), ncInt, ncells_plus_1_dim);

	return boost::bind(&Grid_netcdf_write, this, &nc, vname);
}

/** @param fname Name of file to load from (eg, an overlap matrix file)
@param vname Eg: "grid1" or "grid2" */
void Grid::read_from_netcdf(
NcFile &nc,
std::string const &vname)
{
	clear();

	// ---------- Read the Basic Info
	NcVar *info_var = nc.get_var((vname + ".info").c_str());
		name = std::string(info_var->get_att("name")->as_string(0));
		stype = std::string(info_var->get_att("type")->as_string(0));
		_ncells_full = info_var->get_att("cells.num_full")->as_int(0);
		_nvertices_full = info_var->get_att("vertices.num_full")->as_int(0);

	// ---------- Read the Vertices
	// Basic Info
	std::vector<int> vertices_index(read_int_vector(nc, vname + ".vertices.index"));

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

	// ---------- Read the Cells
	std::vector<int> cells_index(read_int_vector(nc, vname + ".cells.index"));
	std::vector<int> cells_i(read_int_vector(nc, vname + ".cells.i"));
	std::vector<int> cells_j(read_int_vector(nc, vname + ".cells.j"));
	std::vector<int> cells_k(read_int_vector(nc, vname + ".cells.k"));
	std::vector<double> cells_native_area(read_double_vector(nc, vname + ".cells.native_area"));
	std::vector<double> cells_proj_area(read_double_vector(nc, vname + ".cells.proj_area"));

	std::vector<int> vrefs(read_int_vector(nc, vname + ".cells.vertex_refs"));
	std::vector<int> vrefs_start(read_int_vector(nc, vname + ".cells.vertex_refs_start"));

	// Assemble into Cells
	for (size_t i=0; i < cells_index.size(); ++i) {
		int index = cells_index[i];

		Cell cell;
		cell.index = cells_index[i];
		cell.i = cells_i[i];
		cell.j = cells_j[i];
		cell.k = cells_k[i];
		cell.native_area = cells_native_area[i];
		cell.proj_area = cells_proj_area[i];

		// Add the vertices
		cell.reserve(vrefs_start[i+1] - vrefs_start[i]);
		for (int j = vrefs_start[i]; j < vrefs_start[i+1]; ++j)
			cell.add_vertex(get_vertex(vrefs[j]));

		// Add thecell to the grid
		add_cell(std::move(cell));
	}
}

// ---------------------------------------------------
/** @param fname Name of file to load from (eg, an overlap matrix file)
@param vname Eg: "grid1" or "grid2" */
std::unique_ptr<Grid> Grid::netcdf_read(NcFile &nc, std::string const &vname)
{
	auto info_var = nc.get_var((vname + ".info").c_str());
	std::string stype(info_var->get_att("type")->as_string(0));

	std::unique_ptr<Grid> grid;
	if (stype == "xy") {
		grid.reset(new Grid_XY(name, index_base, max_index));
	} else if (stype == "lonlat") {
		grid.reset(new Grid_LatLon(name, index_base, max_index));
	} else {
		grid.reset(new Grid("generic", name, index_base, max_index));
	}


	grid->read_from_netcdf(nc, vname);
	return grid;
}




}	// namespace
