#include <unordered_map>
#include <glint2/exchange_grid.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/bounding_box.h>

namespace glint2 {

/** GLINT2-specific instatiations of CGAL templates. */
namespace gc {
	typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
	typedef Kernel::Point_2                                   Point_2;
	typedef Kernel::Iso_rectangle_2                           Iso_rectangle_2;
	typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
	typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;
}


// Some Supporting Classes
// =======================================================================
static struct OCell {
	Cell *cell;

	/** The polygon representing the grid cell (on the map).
	Vertices in this polygon are always counter-clockwise and have positive area. */
	gc::Polygon_2 poly;
	
	/** Bounding box of the polygon, used for search/overlap algorithms. */
	gc::Iso_rectangle_2 bounding_box;

	OCell(Cell *_cell);
};

OCell::OCell(Cell *_cell) : cell(_cell)
{
	// Copy the vertices
	for (auto vertex = cell->begin(); vertex != cell->end(); +=vertex)
		poly.push_back(gc::Point(vertex->x, vertex->y));

	// Compute the bounding box
	bounding_box = CGAL::bounding_box(poly.vertices_begin(), poly.vertices_end());
}

// =======================================================================

static struct OGrid {
	Grid *grid;

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
	typedef ::RTree<OCell const *, double, 2, double> RTree;
	std::unique_ptr<Grid::RTree> rtree;

	// -------------------------------------------

	OGrid(Grid *_grid) : grid(_grid);

	void realize_rtree();
};		// struct OGrid

OGrid::OGrid(Grid *_grid) : grid(_grid) {

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
		OCell ocell(&*cell);
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
	for (auto ii1=this->cells().begin(); ii1 != this->cells().end(); ++ii1) {
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
		rtree->Insert(min, max, &gc);
	}
}

// =======================================================================
// A general template to compute CGAL Polygon Overlaps

/**
Computes the overlap area of two linear simple polygons. Uses BSO and then add the area of the polygon and substract
the area of its holes.
@param P The first polygon.
@param Q The second polygon.
@return The area of the overlap between P and Q.
@see acg.cs.tau.ac.il/courses/workshop/spring-2008/useful-routines.h
*/
template <class Kernel, class Container>
CGAL::Polygon_2<Kernel, Container>
poly_overlap(const CGAL::Polygon_2<Kernel, Container> &P, 
	const CGAL::Polygon_2<Kernel, Container> &Q)
{
	CGAL_precondition(P.is_simple());
	CGAL_precondition(Q.is_simple());

	// typedef typename CGAL::Polygon_2<Kernel, Container>::FT FT;
	typedef CGAL::Polygon_with_holes_2<Kernel, Container> Polygon_with_holes_2;
	// typedef std::list<Polygon_with_holes_2> Pol_list;

	std::list<Polygon_with_holes_2> overlap;
	CGAL::intersection(P, Q, std::back_inserter(overlap));

	// Look at all the polygons returned
	auto ii0 = overlap.begin();
	if (ii0 == overlap.end()) {
		// Empty list
		return Polygon_2();
	}

	auto ii1 = ii0 + 1;
	if (ii1 != overlap.end()) {
		// More than one Polygon here --- complain
		fprintf(stderr, "ERROR: Overlap expects only simple polygons.  If this is not a bug, you must upgrade your non-CGAL overlap data structures.");
		throw std::err();
	}

	if (ii0->holes_begin() != ii0->holes_end) {
		// Polygon has holes --- again, our format doesn't support it.
		fprintf(stderr, "ERROR: Overlap expects only simple polygons.  If this is not a bug, you must upgrade your non-CGAL overlap data structures.");
		throw std::err();
	}

	// We have a nice, clean overlap.  Return it!
	return ii0->outer_boundary;
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
	OCell **ocell1p, OCell *ocell2)
{
	// Enable using same boost::function callback for many values of grid1
	OCell *ocell1 = *ocell1p;

	// Compute the overlap polygon (CGAL)
	auto expoly(poly_overlap(ocell1->poly, ocell2->poly));
	if (expoly.size() == 0) return true;

	// Convert it to a glint2::Cell
	Cell excell;	// Exchange Cell
	excell.i = ocell1->cell->index;
	excell.j = ocell2->cell->index;
	excell.index = excell.i * grid2_full_nvertices + excell.j;	// guarantee unique

	// Add the vertices of the polygon outline
	for (auto vertex = expoly.begin(); vertex != expoly.end(); ++vertex) {
		double x = CGAL::to_double(vertex->x())
		double y = CGAL::to_double(vertex->y())
		exvcache->add_vertex(excell, x, y);
	}

	// Add it to the grid
	exvcache->grid->add_cell(excell);
}
// --------------------------------------------------------------------

/** @param grid2 Put in an RTree */
std::unique_ptr<Grid> compute_exchange_grid(Grid &grid1, Grid &grid2)
{
	/** Initialize the new grid */
	std::unique_ptr<Grid> exgrid(New Grid("exchange"));
	exgrid->name = grid1.name + '-' + grid2.name;
	exgrid->ncells_full = grid1.ncells_full * grid2.ncells_full;
	exgrid->nvertices_full = -1;	// Not specified
	VertexCache exvcache(VertexCache(&*exgrid));

	OGrid ogrid1(grid1);
	OGrid ogrid2(grid2);
	ogrid2.realize_rtree();

	OCell const *ocell1;
	auto callback(boost::bind(&overlap_callback, &exvcache,
		grid2->full_nvertices, &ocell1, _1));

	for (auto ii1 = ogrid1.ocells.begin(); ii1 != ogrid1.ocells.end(); ++ii1) {
		ocell1 = &ii1->second;		// Set parameter for the callback

		double min[2];
		double max[2];

		min[0] = CGAL::to_double(ocell1->bounding_box.xmin());
		min[1] = CGAL::to_double(ocell1->bounding_box.ymin());
		max[0] = CGAL::to_double(ocell1->bounding_box.xmax());
		max[1] = CGAL::to_double(ocell1->bounding_box.ymax());

		int nfound = ogrid2.rtree->Search(min, max, callback);

		// Logging
		++i;
		if (i % 1000 == 0) {
			printf("Processed %d of %d from grid1, total overlaps = %d\n",
				i+1, ogrid1.cells.size(), exgrid->ncells_realized());
		}
	}

	return exgrid;	
}

};	// namespace glint2
