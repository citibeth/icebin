#pragma once

#include <vector>
#include <unordered_map>

namespace glint2 {

struct Vertex : public Vertex {
	int index;
	double x, double y;
	Vertex(double _x, double _y, _index=-1) : index(_index), x(_x), y(_y) {}

	bool operator<(Cell const &rhs) const
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

	double native_area() { return _native_area; }
	double proj_area() {	// Lazy evaluation here
		if (_proj_area >= 0) { return _proj_area; }

		_proj_area = area_of_polygon(*this);
		return _proj_area;
	}



	size_t size() { return _vertices.size(); }
	typedef DerefIterator<std::vector<Vertex *>> vertex_iterator;
	vertex_iterator begin() { return vertex_iterator(_vertices.begin()); }
	vertex_iterator end() { return vertex_iterator(_vertices.end()); }

	void reserve(size_t n) { _vertices.reserve(n); }
	void add_vertex(Vertex *vertex) { _vertices.push_back(vertex); }

	void clear() {
		_vertices::clear();
		index = -1;
		i = -1;
		j = -1;
		k = -1;
		_native_area = -1;
		_proj_area = -1;
	}

	Cell() { clear(); }
};
// ----------------------------------------------------
class Grid {
	Dict<int, Vertex> _vertices;
	Dict<int, Cell> _cells;

public:
	std::string stype;
	std::string name;

	long ncells_full;		// Maximum possible index (-1)
	long nvertices_full;	// Maximum possible index (-1)

	Grid(std::string const &_stype) : stype(_stype) {}

	// ========= Cell Collection Operators
	// typedef DerefIterator<Dict<int, Cell>> cell_iterator;
	typedef Dict<int, Vertex>::ValIterator cell_iterator;
	iterator cells_begin() { return cell_iterator(_cells.begin()); }
	iterator cells_end() { return cell_iterator(_cells.end()); }

	Cell *get_cell(int index) { return _cells[index]; }
	long ncells_realized() { return _cells.size(); }

	Cell *add_cell(Cell &&cell) {
		// If we never specify our indices, things will "just work"
		if (cell.index == -1) cell.index = _cells.size();


		auto ret = _cells.insert(cell.index, std::move(cell));
		Cell *valp = ret.first;
		bool inserted = ret.second;

		if (!inserted) {		// Key already existed
			fprintf(stderr, "Error adding repeat cell index=%d.  "
				"Cells must have unique indices.", cell.index);
			throw std::exception;
		}
		return valp;
	}

	// ========= Vertex Collection Operators
	typedef Dict<int, Cell>::ValIterator vertex_iterator;
	vertex_iterator vertices_begin() { return vertex_iterator(_vertices.begin()); }
	vertex_iterator vertices_end() { return vertex_iterator(_vertices.end()); }

	Vertex *get_vertex(int index) { return _vertices[index]; }

	long nvertices_realized() { return _vertices.size(); }

	Vertex *add_vertex(Vertex &&vertex) {
		// If we never specify our indices, things will "just work"
		if (vertex.index == -1) vertex.index = _vertices.size();

		auto ret = _vertices.insert(vertex.index, std::move(vertex));
		Vertex *valp = ret.first;
		bool inserted = ret.second;

		if (!inserted) {		// Key already existed
			fprintf(stderr, "Error adding repeat vertex index=%d.  "
				"Vertices must have unique indices.", vertex.index);
			throw std::exception;
		}
		return valp;
	}

	// ========================================


	virtual boost::function<void()> netcdf_define(NcFile &nc, std::string const &vname) const;

	std::unique_ptr<Grid_LonLat> read_from_netcdf(NcFile &nc, std::string const &vname) const;

	std::unique_ptr<Grid> netcdf_read(NcFile &nc, std::string const &vname);

};


}	// namespaces
