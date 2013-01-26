#include <unordered_map>
#include <boost/bind.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/bounding_box.h>

#include <giss/RTree.hpp>
#include <giss/Proj2.hpp>

#include <glint2/ExchangeGrid.hpp>
#include <glint2/gridutil.hpp>

#include <glint2/cgal.hpp>

namespace glint2 {

// Some Supporting Classes
// =======================================================================
struct OCell {
	Cell const *cell;

	/** The polygon representing the grid cell (on the map).
	Vertices in this polygon are always counter-clockwise and have positive area. */
	gc::Polygon_2 poly;
	
	/** Bounding box of the polygon, used for search/overlap algorithms. */
	gc::Iso_rectangle_2 bounding_box;

	OCell(Cell const *_cell, giss::Proj2 const &proj);
};

OCell::OCell(Cell const *_cell, giss::Proj2 const &proj) : cell(_cell)
{
	// Copy the vertices
	for (auto vertex = cell->begin(); vertex != cell->end(); ++vertex) {
		double x, y;
		proj.transform(vertex->x, vertex->y, x, y);
		poly.push_back(gc::Point_2(x, y));
	}

	// Compute the bounding box
	bounding_box = CGAL::bounding_box(poly.vertices_begin(), poly.vertices_end());
}

// =======================================================================

struct OGrid {
	Grid const *grid;

	/** CGAL polygon for each grid cell */
	std::unordered_map<int, OCell> ocells;

	/** Lazily computes an overall bounding box for all realized grid cells.
	The bounding box may or may not be rectangular.
	If the grid has changed since the last time this method was
	called, the bounding box is recomputed.
	@return A polygon containing all realized grid cells. */
	gc::Polygon_2 bounding_box;

	/** RTree index into the Cells of this Grid.
	Allows for efficient searching over limited areas. */
	typedef giss::RTree<OCell const *, double, 2, double> RTree;
	std::unique_ptr<RTree> rtree;

	// -------------------------------------------

	OGrid(Grid const *_grid, giss::Proj2 const &proj);

	void realize_rtree();
};		// struct OGrid

OGrid::OGrid(Grid const *_grid, giss::Proj2 const &proj) : grid(_grid) {

	// Compute bounding box too
	// Be lazy, base bounding box on minimum and maximum values in points
	// (instead of computing the convex hull)
	gc::Kernel::FT minx(1e100);
	gc::Kernel::FT maxx(-1e100);
	gc::Kernel::FT miny(1e100);
	gc::Kernel::FT maxy(-1e100);
	// bounding_box.clear();

	// Compute Simple Bounding Box for overall grid
	for (auto cell = grid->cells_begin(); cell != grid->cells_end(); ++cell) {
		// Convert and copy to the OGrid data structure
		OCell ocell(&*cell, proj);
		ocells.insert(std::make_pair(cell->index, ocell));

		for (auto vertex = ocell.poly.vertices_begin(); vertex != ocell.poly.vertices_end(); ++vertex) {
			minx = std::min(minx, vertex->x());
			maxx = std::max(maxx, vertex->x());
			miny = std::min(miny, vertex->y());
			maxy = std::max(maxy, vertex->y());
		}
	}

	// Store it away
	bounding_box.push_back(gc::Point_2(minx, miny));
	bounding_box.push_back(gc::Point_2(maxx, miny));
	bounding_box.push_back(gc::Point_2(maxx, maxy));
	bounding_box.push_back(gc::Point_2(minx, maxy));
}

void OGrid::realize_rtree() {
	rtree.reset(new OGrid::RTree);

	double min[2];
	double max[2];
	for (auto ii1=ocells.begin(); ii1 != ocells.end(); ++ii1) {
		OCell &ocell(ii1->second);

		min[0] = CGAL::to_double(ocell.bounding_box.xmin());
		min[1] = CGAL::to_double(ocell.bounding_box.ymin());
		max[0] = CGAL::to_double(ocell.bounding_box.xmax());
		max[1] = CGAL::to_double(ocell.bounding_box.ymax());

		//fprintf(stderr, "Adding bounding box: (%f %f)  (%f %f)\n", min[0], min[1], max[0], max[1]);

		// Deal with floating point...
		const double eps = 1e-7;
		double epsilon_x = eps * std::abs(max[0] - min[0]);
		double epsilon_y = eps * std::abs(max[1] - min[1]);
		min[0] -= epsilon_x;
		min[1] -= epsilon_y;
		max[0] += epsilon_x;
		max[1] += epsilon_y;

		//std::cout << ocell.poly << std::endl;
		//std::cout << ocell.bounding_box << std::endl;
		//printf("(%g,%g) -> (%g,%g)\n", min[0], min[1], max[0], max[1]);
		rtree->Insert(min, max, &ocell);
	}
}


// =======================================================================
// The main exchange grid computation

/**
@param exgrid The Exchange Grid we're creating.  Even for L1 grids, we
	don't need to positively associate vertices in exgrid with
	vertices in grid1 or grid2.  No more than a "best effort" is
	needed to eliminate duplicate vertices.
@return Always returns true (tells RTree search algorithm to keep going) */
static bool overlap_callback(VertexCache *exvcache, long grid2_full_nvertices,
	OCell const **ocell1p, OCell const *ocell2)
{
	// Enable using same boost::function callback for many values of grid1
	OCell const *ocell1 = *ocell1p;

	// Compute the overlap polygon (CGAL)
	auto expoly(poly_overlap(ocell1->poly, ocell2->poly));
	if (expoly.size() == 0) return true;

	// Convert it to a glint2::Cell
	Cell excell;	// Exchange Cell
	excell.i = ocell1->cell->index;
	excell.j = ocell2->cell->index;
	excell.index = excell.i * grid2_full_nvertices + excell.j;	// guarantee unique

	// Add the vertices of the polygon outline
	for (auto vertex = expoly.vertices_begin(); vertex != expoly.vertices_end(); ++vertex) {
		double x = CGAL::to_double(vertex->x());
		double y = CGAL::to_double(vertex->y());
		exvcache->add_vertex(excell, x, y);
	}

	// Compute its area
	excell.area = area_of_polygon(excell);

	// Add it to the grid
	exvcache->grid->add_cell(std::move(excell));

	return true;
}
// --------------------------------------------------------------------

/** @param grid2 Put in an RTree */
//std::unique_ptr<Grid> compute_exchange_grid
ExchangeGrid::ExchangeGrid(Grid const &grid1, Grid const &grid2, std::string const &_sproj)
: Grid("exchange")
{
	scoord = "xy";
	ExchangeGrid *exgrid = this;

printf("ExchangeGrid 1\n");

	/** Suss out projections */
printf("grid1.scoord = %s, grid2.scoord = %s\n", grid1.scoord.c_str(), grid2.scoord.c_str());
printf("grid1.sproj = %s, grid2.sproj = %s\n", grid1.sproj.c_str(), grid2.sproj.c_str());

	if (grid1.scoord == "xy") {
		if (grid2.scoord == "xy") {
			// No projections needed
			if (grid1.sproj != grid2.sproj) {
				fprintf(stderr, "Two XY grids must have the same projection\n");
				throw std::exception();
			}
		} else {
			// grid1=xy, grid2=ll: Project from grid 2 to grid1's xy
			proj2.init(grid1.sproj, giss::Proj2::Direction::LL2XY);
			sproj = (_sproj == "" ? std::string(grid1.sproj.c_str()) : _sproj);
		}
	} else {
		if (grid2.scoord == "xy") {
			// grid1=ll, grid2=xy: Project from grid 1 to grid2's xy
			proj1.init(grid2.sproj, giss::Proj2::Direction::LL2XY);
			sproj = (_sproj == "" ? std::string(grid2.sproj.c_str()) : _sproj);
		} else {
			// Both in Lat/Lon: Project them both to XY for overlap computation
			// BUT... since we don't have a projection, we must throw
			// up our hands!
			fprintf(stderr, "Program isn't currently equipped to overlap two grids on the sphere.\n");
			throw std::exception();
		}
	}

printf("ExchangeGrid 2\n");

	/** Initialize the new grid */
	exgrid->name = grid1.name + '-' + grid2.name;
	exgrid->grid1_ncells_full = grid1.ncells_full;
	exgrid->grid2_ncells_full = grid2.ncells_full;
	exgrid->ncells_full = grid1.ncells_full * grid2.ncells_full;
	exgrid->nvertices_full = -1;	// Not specified
	VertexCache exvcache(exgrid);

	OGrid ogrid1(&grid1, proj1);
	OGrid ogrid2(&grid2, proj2);
	ogrid2.realize_rtree();

	OCell const *ocell1;
	auto callback(boost::bind(&overlap_callback, &exvcache,
		grid2.nvertices_full, &ocell1, _1));

	int nprocessed=0;
	for (auto ii1 = ogrid1.ocells.begin(); ii1 != ogrid1.ocells.end(); ++ii1) {
		ocell1 = &ii1->second;		// Set parameter for the callback

		double min[2];
		double max[2];

printf("grid1[%d]: x in (%f - %f), y in (%f - %f)\n", ocell1->cell->index, min[0], max[0], min[1], max[1]);

		min[0] = CGAL::to_double(ocell1->bounding_box.xmin());
		min[1] = CGAL::to_double(ocell1->bounding_box.ymin());
		max[0] = CGAL::to_double(ocell1->bounding_box.xmax());
		max[1] = CGAL::to_double(ocell1->bounding_box.ymax());

		int nfound = ogrid2.rtree->Search(min, max, callback);

		// Logging
		++nprocessed;
		if (nprocessed % 1000 == 0) {
			printf("Processed %d of %d from grid1, total overlaps = %d\n",
				nprocessed+1, ogrid1.ocells.size(), exgrid->ncells_realized());
		}
	}
}

// ---------------------------------------------------------------
boost::function<void ()> ExchangeGrid::netcdf_define(NcFile &nc, std::string const &vname) const
{
	auto parent = Grid::netcdf_define(nc, vname);

	NcVar *info_var = nc.get_var((vname + ".info").c_str());
	info_var->add_att("grid1.ncells_full", grid1_ncells_full);
	info_var->add_att("grid2.ncells_full", grid2_ncells_full);
	proj1.netcdf_define(nc, info_var, "proj1");
	proj2.netcdf_define(nc, info_var, "proj2");

	return parent;
}

void ExchangeGrid::read_from_netcdf(NcFile &nc, std::string const &vname)
{
	Grid::read_from_netcdf(nc, vname);

	NcVar *info_var = nc.get_var((vname + ".info").c_str());
	grid1_ncells_full = info_var->get_att("grid1.ncells_full")->as_int(0);
	grid2_ncells_full = info_var->get_att("grid2.ncells_full")->as_int(0);
	proj1.read_from_netcdf(nc, info_var, "proj1");
	proj2.read_from_netcdf(nc, info_var, "proj2");
}


};	// namespace glint2
