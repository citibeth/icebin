#include <cmath>
#include <boost/bind.hpp>
#include <giss/ncutil.hpp>
#include <glint2/Grid_XY.hpp>
#include <glint2/gridutil.hpp>

namespace glint2 {

void set_xy_boundaries(Grid_XY &grid,
	double x0, double x1, double dx,
	double y0, double y1, double dy)
{
	// Set up x coordinates
	grid.xb.clear();
	int nx = (int)(.5 + (x1 - x0) / dx);	// Round to nearest integer
	double nx_inv = 1.0 / (double)nx;
	for (int i=0; i<nx; ++i) {
		double x = x0 + (x1-x0) * (double)i * nx_inv;
		grid.xb.push_back(x);
	}

	// Set up y coordinates
	grid.yb.clear();
	int ny = (int)(.5 + (y1 - y0) / dy);	// Round to nearest integer
	double ny_inv = 1.0 / (double)ny;
	for (int i=0; i<ny; ++i) {
		double y = y0 + (y1-y0) * (double)i * ny_inv;
		grid.yb.push_back(y);
	}

}

void set_xy_centers(Grid_XY &grid,
	double x0, double x1, double dx,
	double y0, double y1, double dy)
{
	set_xy_boundaries(grid,
		x0-.5*dx, x1+.5*dx, dx,
		y0-.5*dy, y1+.5*dy, dy);
}


void Grid_XY::realize_grid(
boost::function<bool(Cell const &)> const &euclidian_clip)
{
	ncells_full = nx() * ny();
	nvertices_full = xb.size() * yb.size();

	// Set up the main grid
	VertexCache vcache(this);
	int index = 1;
	for (int iy = 0; iy < yb.size()-1; ++iy) {
		double y0 = yb[iy];
		double y1 = yb[iy+1];
//		grid->y_centers.push_back(.5*(y0+y1));

		for (int ix = 0; ix < xb.size()-1; ++ix, ++index) {
			double x0 = xb[ix];
			double x1 = xb[ix+1];
//			grid->x_centers.push_back(.5*(x0+x1));

			Cell cell;
			vcache.add_vertex(cell, x0, y0);
			vcache.add_vertex(cell, x1, y0);
			vcache.add_vertex(cell, x1, y1);
			vcache.add_vertex(cell, x0, y1);

			// Don't include things outside our clipping region
			if (!euclidian_clip(cell)) continue;

			cell.index = index;
			cell.i = ix;
			cell.j = iy;
			cell._native_area = area_of_polygon(cell);
			cell._proj_area = cell._native_area;

			add_cell(std::move(cell));
		}
	}
}

static void Grid_XY_netcdf_write(
	boost::function<void()> const &parent,
	NcFile *nc, Grid_XY const *grid, std::string const &vname)
{
	parent();

	NcVar *xbVar = nc->get_var((vname + ".x_boundaries").c_str());
	NcVar *ybVar = nc->get_var((vname + ".y_boundaries").c_str());

	xbVar->put(&grid->xb[0], grid->xb.size());
	ybVar->put(&grid->yb[0], grid->yb.size());
}

boost::function<void ()> Grid_XY::netcdf_define(NcFile &nc, std::string const &vname) const
{
	auto parent = Grid::netcdf_define(nc, vname);

	NcDim *xbDim = nc.add_dim((vname + ".x_boundaries.length").c_str(),
		this->xb.size());
	NcVar *xbVar = nc.add_var((vname + ".x_boundaries").c_str(),
		ncDouble, xbDim);
	NcDim *ybDim = nc.add_dim((vname + ".y_boundaries.length").c_str(),
		this->yb.size());
	NcVar *ybVar = nc.add_var((vname + ".y_boundaries").c_str(),
		ncDouble, ybDim);

	return boost::bind(&Grid_XY_netcdf_write, parent, &nc, this, vname);
}

void Grid_XY::read_from_netcdf(NcFile &nc, std::string const &vname)
{
	Grid::read_from_netcdf(nc, vname);

	xb = giss::read_double_vector(nc, vname + ".x_boundaries");
	yb = giss::read_double_vector(nc, vname + ".y_boundaries");
}


}
