#pragma once

#include <vector>
#include <unordered_map>
#include <netcdfcpp.h>
#include <boost/function.hpp>
#include <giss/Dict.hpp>
#include <giss/Proj2.hpp>
#include <giss/enum.hpp>

namespace glint2 {

class Cell;

/** See Surveyor's Formula: http://www.maa.org/pubs/Calc_articles/ma063.pdf*/
extern double area_of_polygon(Cell const &cell);
/** @param proj[2] */
extern double area_of_proj_polygon(Cell const &cell, giss::Proj2 const &proj);

// --------------------------------------------------

struct Vertex {
	long index;
	double const x;
	double const y;

	Vertex(double _x, double _y, int _index=-1) :
		x(_x), y(_y), index(_index) {}

	bool operator<(Vertex const &rhs) const
		{ return index < rhs.index; }
};

// ----------------------------------------------------
/** Iterate through with:
<pre>Cell cell;
for (auto ii=cell.begin(); ii != cell.end(); ++ii) {
	print ("Vertex %d\n", ii->index);
}</pre>
*/
class Cell {
	std::vector<Vertex *> _vertices;

public:

	/** For L0 formulations (constant value per grid cell):
	Index of this grid cell in dense arrays (base=0) */
	long index;

	bool operator<(Cell const &rhs) const
		{ return index < rhs.index; }

	/** Optional stuff.
	For exchange grid: i, j tell the source coordinates (0-based)
	For grids with 2-D indexing, tells the i and j index of the cell (0-based). */
	int i, j, k;

	/** Area of this grid cell, whether it is a Cartesian or lat/lon cell. */
	double area;

	size_t size() const { return _vertices.size(); }
	typedef giss::DerefIterator<std::vector<Vertex *>::const_iterator> vertex_iterator;
	vertex_iterator begin() const { return vertex_iterator(_vertices.begin()); }
	vertex_iterator end() const { return vertex_iterator(_vertices.end()); }

	void reserve(size_t n) { _vertices.reserve(n); }
	void add_vertex(Vertex *vertex) { _vertices.push_back(vertex); }

	void clear() {
		_vertices.clear();
		index = -1;
		i = -1;
		j = -1;
		k = -1;
		area = -1;
//		_native_area = -1;
//		_proj_area = -1;
	}

	Cell() { clear(); }
};		// class Cell
// ----------------------------------------------------
class Grid {
	giss::HashDict<int, Vertex> _vertices;
	giss::HashDict<int, Cell> _cells;

public:
	// Corresponds to classes
	BOOST_ENUM_VALUES( Type, int,
		(GENERIC)		(0)		// Just use the Grid base class
		(XY)			(1)		// Rectilinear X/Y grid
		(LONLAT)		(2)		// Global lat-lon grid (maybe with polar caps)
		(EXCHANGE)		(3)		// Exchange grid, from overlap of two other grids
//		(CUBESPHERE)	(4)		// Global Cubed Sphere grid
//		(MESH)			(5)		// Arbitrary mesh (could be global or on plane)
	)

	BOOST_ENUM_VALUES( Coordinates, int,
		(XY)			(0)		// Vertices in x/y coordinates on a plane
		(LONLAT)		(1)		// Vertices in lon/lat coordinates on a sphere
	)

	BOOST_ENUM_VALUES( Parameterization, int,
		(L0)			(0)		// Constant value in each grid cell
		(L1)			(1)		// Value specified at each vertex, slope inbetween
	)

	Type type;
	Coordinates coordinates;
	Parameterization parameterization;

	std::string name;

protected :
	long _ncells_full;		// Maximum possible index (-1)
	long _nvertices_full;	// Maximum possible index (-1)

	// These are kept in line, with add_cell() and add_vertex()
	long _max_realized_cell_index;		// Maximum index of realized cells
	long _max_realized_vertex_index;
public:

	long ncells_full() const
		{ return _ncells_full >= 0 ? _ncells_full : _max_realized_cell_index; }
	long nvertices_full() const
		{ return _nvertices_full >= 0 ? _nvertices_full : _max_realized_vertex_index; }

	/** Projection, if any, used to transform this grid from its
	"native" form to the one it's currently in.
	Does transform(proj[0], proj[1]).
	sproj[0] == sproj[1] == "" if no projection used. */
//	std::string srpoj[2];

	/** If scoord == "xy": The projection that relates x,y coordinates
	here to a specific point on the globe (as a Proj.4 String). */
	std::string sproj;

	Grid(Type type);
	virtual ~Grid() {}

	/** @return ncells_full (for L0) or nvertices_full (for L1) */
	long ndata() const;

	virtual void clear();

	/** For regularly spaced grids: converts 2D (or 3D in some cases) indexing to 1-D.
	The index will be a Cell index for L0 grids, or a Vertex index for L1 */
	virtual int ij_to_index(int i, int j) const { return -1; }

	virtual void index_to_ij(long index, int &i, int &j) const { }

	// ========= Cell Collection Operators
	auto cells_begin() -> decltype(_cells.begin()) { return _cells.begin(); }
	auto cells_end() -> decltype(_cells.end()) { return _cells.end(); }

	auto cells_begin() const -> decltype(_cells.begin()) { return _cells.begin(); }
	auto cells_end() const -> decltype(_cells.end()) { return _cells.end(); }

	Cell *get_cell(long index) { return _cells[index]; }
	long ncells_realized() const { return _cells.size(); }

	Cell *add_cell(Cell &&cell);

	void cells_erase(giss::HashDict<int, Cell>::iterator &ii)
		{ _cells.erase(ii); }

	// ========= Vertex Collection Operators

	auto vertices_begin() -> decltype(_vertices.begin()) { return _vertices.begin(); }
	auto vertices_end() -> decltype(_vertices.end()) { return _vertices.end(); }

	auto vertices_begin() const -> decltype(_vertices.begin())
		{ return _vertices.begin(); }
	auto vertices_end() const -> decltype(_vertices.end())
		{ return _vertices.end(); }

	Vertex *get_vertex(long index) { return _vertices[index]; }

	long nvertices_realized() const { return _vertices.size(); }

	Vertex *add_vertex(Vertex &&vertex);

	/** Sort vertices by X/Y coordinates,
	and then re-number them into sorted order. */
	void sort_renumber_vertices();

	// ========================================

	void get_ll_to_xy(giss::Proj2 &proj, std::string const &sproj) const;
	std::vector<double> get_proj_areas(std::string const &sproj) const;
	std::vector<double> get_native_areas() const;



protected:
	void netcdf_write(NcFile *nc, std::string const &vname) const;

public:
	virtual boost::function<void()> netcdf_define(NcFile &nc, std::string const &vname) const;
	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);

	void to_netcdf(std::string const &fname);


	/** Remove cells and vertices not relevant to us --- for example, not in our MPI domain.
	This will be done AFTER we read it in.  It's an optimization. */
	void filter_cells(boost::function<bool (int)> const &include_cell);

};

std::unique_ptr<Grid> read_grid(NcFile &nc, std::string const &vname);



}	// namespaces
