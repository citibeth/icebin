#pragma once

#include <vector>
#include <unordered_map>
#include <netcdfcpp.h>
#include <boost/function.hpp>
#include <giss/Dict.hpp>

namespace glint2 {

class Cell;

/** See Surveyor's Formula: http://www.maa.org/pubs/Calc_articles/ma063.pdf */
extern double area_of_polygon(Cell const &cell);

// --------------------------------------------------

struct Vertex {
	int index;
	double x;
	double y;
	Vertex(double _x, double _y, int _index=-1) : index(_index), x(_x), y(_y) {}

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
	int index;

	bool operator<(Cell const &rhs) const
		{ return index < rhs.index; }

	/** Optional stuff.
	For exchange grid: i, j tell the source coordinates (0-based)
	For grids with 2-D indexing, tells the i and j index of the cell (0-based). */
	int i, j, k;

	/** Area of this grid cel in its native coordinate system (if it's
	been projected) */
	double _native_area;
	double _proj_area;

	double native_area() const { return _native_area; }
	double proj_area() const {	// Lazy evaluation here
		if (_proj_area >= 0) { return _proj_area; }

		// We allow a const_cast here becuase _proj_area is used for
		// caching values ONLY.
		const_cast<Cell *>(this)->_proj_area = area_of_polygon(*this);
		return _proj_area;
	}



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
		_native_area = -1;
		_proj_area = -1;
	}

	Cell() { clear(); }
};		// class Cell
// ----------------------------------------------------
class Grid {
	giss::Dict<int, Vertex> _vertices;
	giss::Dict<int, Cell> _cells;

public:
	std::string stype;
	std::string name;

	long ncells_full;		// Maximum possible index (-1)
	long nvertices_full;	// Maximum possible index (-1)

	Grid(std::string const &_stype);

	virtual void clear();


	/** For regularly spaced grids: converts 2D (or 3D in some cases) indexing to 1-D.
	The index will be a Cell index for L0 grids, or a Vertex index for L1 */
	virtual int ij_to_index(int i, int j) const { return -1; }

	virtual void index_to_ij(int index, int &i, int &j) const { }

	// ========= Cell Collection Operators
	// typedef DerefIterator<giss::Dict<int, Cell>> cell_iterator;
	auto cells_begin() -> decltype(_cells.begin()) { return _cells.begin(); }
	auto cells_end() -> decltype(_cells.end()) { return _cells.end(); }

	auto cells_begin() const -> decltype(_cells.begin()) { return _cells.begin(); }
	auto cells_end() const -> decltype(_cells.end()) { return _cells.end(); }



#if 0
	typedef giss::Dict<int, Vertex>::ValIterator cell_iterator;
	cell_iterator cells_begin() { return cell_iterator(_cells.begin()); }
	cell_iterator cells_end() { return cell_iterator(_cells.end()); }
#endif

	Cell *get_cell(int index) { return _cells[index]; }
	long ncells_realized() const { return _cells.size(); }

	Cell *add_cell(Cell &&cell);

	// ========= Vertex Collection Operators

	auto vertices_begin() -> decltype(_vertices.begin()) { return _vertices.begin(); }
	auto vertices_end() -> decltype(_vertices.end()) { return _vertices.end(); }

	auto vertices_begin() const -> decltype(_vertices.begin())
		{ return _vertices.begin(); }
	auto vertices_end() const -> decltype(_vertices.end())
		{ return _vertices.end(); }

#if 0
	typedef giss::Dict<int, Cell>::ValIterator vertex_iterator;
	vertex_iterator vertices_begin() { return vertex_iterator(_vertices.begin()); }
	vertex_iterator vertices_end() { return vertex_iterator(_vertices.end()); }
#endif

	Vertex *get_vertex(int index) { return _vertices[index]; }

	long nvertices_realized() const { return _vertices.size(); }

	Vertex *add_vertex(Vertex &&vertex);

	// ========================================


protected:
	void netcdf_write(NcFile *nc, std::string const &vname) const;

public:
	virtual boost::function<void()> netcdf_define(NcFile &nc, std::string const &vname) const;
	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);
};

std::unique_ptr<Grid> read_grid(NcFile &nc, std::string const &vname);



}	// namespaces
