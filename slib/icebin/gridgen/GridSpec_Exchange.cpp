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

#include <unordered_map>
#include <functional>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/bounding_box.h>

#include <ibmisc/RTree.hpp>
#include <ibmisc/Proj2.hpp>
#include <ibmisc/netcdf.hpp>

#include <icebin/error.hpp>

#include <icebin/gridgen/cgal.hpp>
#include <icebin/gridgen/GridSpec_Exchange.hpp>
#include <icebin/gridgen/gridutil.hpp>

using namespace ibmisc;
using namespace std::placeholders;  // for _1, _2, _3...

namespace icebin {

// Some Supporting Classes
// =======================================================================
struct OCell {
	Cell const *cell;

	/** The polygon representing the grid cell (on the map).
	Vertices in this polygon are always counter-clockwise and have positive area. */
	gc::Polygon_2 poly;
	
	/** Bounding box of the polygon, used for search/overlap algorithms. */
	gc::Iso_rectangle_2 bounding_box;

	OCell(Cell const *_cell, Proj2 const *proj);
};

OCell::OCell(Cell const *_cell, Proj2 const *proj) : cell(_cell)
{
	// Copy the vertices
	for (auto vertex = cell->begin(); vertex != cell->end(); ++vertex) {
		double x, y;
		if (proj) {
			proj->transform(vertex->x, vertex->y, x, y);
		} else {
			x = vertex->x;
			y = vertex->y;
		}
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
	typedef ibmisc::RTree<OCell const *, double, 2, double> RTree;
	std::unique_ptr<RTree> rtree;

	// -------------------------------------------

	OGrid(Grid const *_grid, Proj2 const *proj);

	void realize_rtree();
};		// struct OGrid

OGrid::OGrid(Grid const *_grid, Proj2 const *proj) : grid(_grid)
{

	// Compute bounding box too
	// Be lazy, base bounding box on minimum and maximum values in points
	// (instead of computing the convex hull)
	gc::Kernel::FT minx(1e100);
	gc::Kernel::FT maxx(-1e100);
	gc::Kernel::FT miny(1e100);
	gc::Kernel::FT maxy(-1e100);
	// bounding_box.clear();

	// Compute Simple Bounding Box for overall grid
	for (auto cell = grid->cells.begin(); cell != grid->cells.end(); ++cell) {
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
	vertices in gridA or gridI.  No more than a "best effort" is
	needed to eliminate duplicate vertices.
@return Always returns true (tells RTree search algorithm to keep going) */
static bool overlap_callback(VertexCache *exvcache, long gridI_ndata,
	OCell const **ocell1p, OCell const *ocell2)
{
	// Enable using same std::function callback for many values of gridA
	OCell const *ocell1 = *ocell1p;

	// Compute the overlap polygon (CGAL)
	auto expoly(poly_overlap(ocell1->poly, ocell2->poly));
	if (expoly.size() == 0) return true;

	// Convert it to a Cell
	Cell excell;	// Exchange Cell
	excell.i = ocell1->cell->index;
	excell.j = ocell2->cell->index;
//	excell.index = excell.i * gridI_ndata + excell.j;	// guarantee unique (but sparse)
	excell.index = -1;		// Get an index assigned (but dense)...

	// Add the vertices of the polygon outline
	for (auto vertex = expoly.vertices_begin(); vertex != expoly.vertices_end(); ++vertex) {
		double x = CGAL::to_double(vertex->x());
		double y = CGAL::to_double(vertex->y());
		exvcache->add_vertex(excell, x, y);
	}

	// Compute its area (we will need this)
	excell.native_area = excell.proj_area(NULL);

	// Add it to the grid
	exvcache->grid->cells.add(std::move(excell));

	return true;
}
// --------------------------------------------------------------------

/** @param gridI Put in an RTree */
void GridSpec_Exchange::make_grid(Grid &exgrid)
{
	exgrid.coordinates = Grid::Coordinates::XY;
	exgrid.parameterization = Grid::Parameterization::L0;	// Why not?

	// Determine compatibility and projections between the two grids
	std::unique_ptr<Proj2> proj1, proj2;
	if (gridA->coordinates == Grid::Coordinates::XY) {
		if (gridI->coordinates == Grid::Coordinates::XY) {
			// No projections needed
			if (gridA->sproj != gridI->sproj) {
				(*icebin_error)(-1, "Two XY grids must have the same projection\n");
			}
		} else {
			// gridA=xy, gridI=ll: Project from grid 2 to gridA's xy
			proj2.reset(new Proj2(gridA->sproj, Proj2::Direction::LL2XY));
			if (sproj == "") sproj = std::string(gridA->sproj.c_str());
		}
	} else {
		if (gridI->coordinates == Grid::Coordinates::XY) {
			// gridA=ll, gridI=xy: Project from grid 1 to gridI's xy
			proj1.reset(new Proj2(gridI->sproj, Proj2::Direction::LL2XY));
			if (sproj == "") sproj = std::string(gridI->sproj.c_str());
		} else {
			// Both in Lat/Lon: Project them both to XY for overlap computation
			// BUT... since we don't have a projection, we must throw
			// up our hands!
			(*icebin_error)(-1, "Program isn't currently equipped to overlap two grids on the sphere.\n");
		}
	}

	/** Initialize the new grid */
	exgrid.name = gridA->name + '-' + gridI->name;
	gridA_cells_nfull = gridA->cells.nfull();
	gridI_cells_nfull = gridI->cells.nfull();
	exgrid.cells._nfull = -1;		// Not specified
	exgrid.vertices._nfull = -1;	// Not specified
	VertexCache exvcache(&exgrid);

	OGrid ogridA(gridA, &*proj1);
	OGrid ogridI(gridI, &*proj2);
	ogridI.realize_rtree();

	OCell const *ocell1;
	auto callback(std::bind(&overlap_callback, &exvcache,
		gridI->ndata(), &ocell1, _1));

	int nprocessed=0;
	for (auto ii1 = ogridA.ocells.begin(); ii1 != ogridA.ocells.end(); ++ii1) {
		ocell1 = &ii1->second;		// Set parameter for the callback

		double min[2];
		double max[2];

//printf("gridA[%d]: x in (%f - %f), y in (%f - %f)\n", ocell1->cell->index, min[0], max[0], min[1], max[1]);

		min[0] = CGAL::to_double(ocell1->bounding_box.xmin());
		min[1] = CGAL::to_double(ocell1->bounding_box.ymin());
		max[0] = CGAL::to_double(ocell1->bounding_box.xmax());
		max[1] = CGAL::to_double(ocell1->bounding_box.ymax());

		int nfound = ogridI.rtree->Search(min, max, callback);

		// Logging
		++nprocessed;
		if (nprocessed % 10 == 0) {
			printf("Processed %d of %d from gridA, total overlaps = %d\n",
				nprocessed+1, ogridA.ocells.size(), exgrid.cells.nrealized());
		}
	}

	exgrid.indexing = Indexing<int, long>({0}, {exgrid.cells.nfull()}, {0});	// No n-D indexing available.

}


};	// namespace glint2
