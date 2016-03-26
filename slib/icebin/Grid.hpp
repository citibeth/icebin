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

#pragma once

#include <vector>
#include <unordered_map>
#include <functional>

#include <ibmisc/enum.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/iter.hpp>
#include <ibmisc/Proj2.hpp>
#include <ibmisc/indexing.hpp>

#include <icebin/error.hpp>
#include <icebin/sparse.hpp>

namespace icebin {

class Cell;

// --------------------------------------------------

struct Vertex {
    long index;
    double const x;
    double const y;

    Vertex() : x(0), y(0), index(-1) {}

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

    /** Area of this grid cell, in its native (unprojected) coordinate
    system. */
    double native_area;

    Cell() : index(-1), i(-1), j(-1), k(-1), native_area(0) {}

    Cell(std::vector<Vertex *> &&vertices) : 
        _vertices(std::move(vertices)),
        index(-1), i(-1), j(-1), k(-1), native_area(0) {}

    Cell(std::vector<Vertex *> &vertices) : 
        _vertices(vertices),
        index(-1), i(-1), j(-1), k(-1), native_area(0) {}

    bool operator<(Cell const &rhs) const
        { return index < rhs.index; }

    /** Optional stuff.
    For exchange grid: i, j tell the source coordinates (0-based)
    For grids with 2-D indexing, tells the i and j index of the cell (0-based). */
    int i, j, k;

    size_t size() const { return _vertices.size(); }

    typedef ibmisc::DerefRandomAccessIter<Vertex, std::vector<Vertex *>::iterator> iterator;
    typedef ibmisc::DerefRandomAccessIter<const Vertex, std::vector<Vertex *>::const_iterator> const_iterator;

    iterator begin(int ix = 0)
        { return iterator(_vertices.begin() + ix); }
    iterator end(int ix = 0)
        { return iterator(_vertices.end() + ix); }
    const_iterator cbegin(int ix = 0) const
        { return const_iterator(_vertices.cbegin() + ix); }
    const_iterator cend(int ix = 0) const
        { return const_iterator(_vertices.cend() + ix); }
    const_iterator begin(int ix = 0) const
        { return const_iterator(_vertices.cbegin() + ix); }
    const_iterator end(int ix = 0) const
        { return const_iterator(_vertices.cend() + ix); }

    void reserve(size_t n) { _vertices.reserve(n); }
    void add_vertex(Vertex *vertex) { _vertices.push_back(vertex); }

    double proj_area(ibmisc::Proj_LL2XY const *proj);   // OPTIONAL

};      // class Cell
// ----------------------------------------------------
class Grid;

class GridSpec {};  // Tagging class

/** Specialized dict-like structure used for cells and vertices in a grid. */
template<class CellT>
class GridMap {
    friend class Grid;
    friend class GridSpec_XY;
    friend class GridSpec_LonLat;
    friend class GridSpec_Exchange;
protected:
    typedef std::unordered_map<long, std::unique_ptr<CellT>> MapT;
    MapT _cells;
    long _nfull = -1;
    long _max_realized_index = -1;

#if 0
    // ----------------
    // Disallow public copying, makin this just a part of the Grid class.
    // http://www.stroustrup.com/C++11FAQ.html#default
    GridMap<CellT>& operator=(const GridMap<CellT>&) = delete;
    GridMap<CellT>(const GridMap<CellT>&) = delete;

    // ----------------
#endif

public:

    typedef ibmisc::DerefSecondIter<long, CellT, typename MapT::iterator> iterator;
    typedef ibmisc::DerefSecondIter<long, const CellT, typename MapT::const_iterator> const_iterator;


    iterator begin()
        { return iterator(_cells.begin()); }
    iterator end()
        { return iterator(_cells.end()); }
    const_iterator cbegin() const
        { return const_iterator(_cells.cbegin()); }
    const_iterator cend() const
        { return const_iterator(_cells.cend()); }
    const_iterator begin() const
        { return const_iterator(_cells.cbegin()); }
    const_iterator end() const
        { return const_iterator(_cells.cend()); }

    iterator erase(iterator const &ii)
        { return iterator(_cells.erase(ii.wrapped)); }

    void clear() { _cells.clear(); }

    CellT *at(long index) { return &*_cells.at(index); }
    CellT const *at(long index) const { return &*_cells.at(index); }
    size_t nrealized() const { return _cells.size(); }
    size_t nfull() const { return _nfull >=0 ? _nfull : _max_realized_index+1; }

    /** Adds a cell and owns it. */
    CellT *add(CellT &&cell);


    /** Adds the item to our collection, then deletes the pointer.
    This is for use with Cython, which doesn't like RValue references. */
    CellT *add_claim(CellT *cell);

private :
    struct CmpPointers {
        bool operator()(CellT const *a, CellT const *b) { return *a < *b; }
    };

public :

    /** @return A sorted vector of (simple pointers to) the values stored in the Dict. */
    std::vector<CellT const *> sorted() const;
};


template<class CellT>
std::vector<CellT const *> GridMap<CellT>::sorted() const
{
    // Make a vector of pointers
    std::vector<CellT const *> ret;
    for (auto ii = begin(); ii != end(); ++ii) ret.push_back(&*ii);

    // Sort the vector      
    std::sort(ret.begin(), ret.end(), CmpPointers());

    return ret;
}   

template<class CellT>
CellT *GridMap<CellT>::add(CellT &&cell)
{
    // If we never specify our indices, things will "just work"
    if (cell.index < 0) cell.index = _cells.size();
    _max_realized_index = std::max(_max_realized_index, cell.index);

    std::unique_ptr<CellT> ptr(new CellT(std::move(cell)));
    auto ret = _cells.insert(std::make_pair(cell.index, std::move(ptr)));
    CellT *valp = ret.first->second.get();
    bool inserted = ret.second;

    if (!inserted) {        // Key already existed
        (*icebin_error)(-1, "Error adding repeat cell/vertex index=%d.  "
            "Cells and Vertices must have unique indices.", cell.index);
    }
    return valp;
}

template<class CellT>
CellT *GridMap<CellT>::add_claim(CellT *cell)
{
    // Make sure we get deleted, even in face of exceptions
    std::unique_ptr<CellT> pcell(cell);
    return add(std::move(*pcell));
}


// ----------------------------------------------------
class Grid {
public:
    GridMap<Vertex> vertices;
    GridMap<Cell> cells;

    // Corresponds to classes
    BOOST_ENUM_VALUES( Type, int,
        (GENERIC)       (0)     // Just use the Grid base class
        (XY)            (1)     // Rectilinear X/Y grid
        (LONLAT)        (2)     // Global lat-lon grid (maybe with polar caps)
        (EXCHANGE)      (3)     // Exchange grid, from overlap of two other grids
//      (CUBESPHERE)    (4)     // Global Cubed Sphere grid
//      (MESH)          (5)     // Arbitrary mesh (could be global or on plane)
    )

    BOOST_ENUM_VALUES( Coordinates, int,
        (XY)            (0)     // Vertices in x/y coordinates on a plane
        (LONLAT)        (1)     // Vertices in lon/lat coordinates on a sphere
    )

    BOOST_ENUM_VALUES( Parameterization, int,
        (L0)            (0)     // Constant value in each grid cell
        (L1)            (1)     // Value specified at each vertex, slope inbetween
    )

    Type type;
    Coordinates coordinates;
    Parameterization parameterization;

    /** Conversion between n-dimensional indexing used natively on the
    grid, and 1-D indexing used in IceBin. */
    ibmisc::Indexing<int, long> indexing;

    std::string name;

protected :
    // These are kept in line, with add_cell() and add_vertex()
    long _max_realized_cell_index;      // Maximum index of realized cells
    long _max_realized_vertex_index;
public:

    /** If scoord == "xy": The projection that relates x,y coordinates
    here to a specific point on the globe (as a Proj.4 String). */
    std::string sproj;

    Grid();
    virtual ~Grid() {}

    /** The size of the vector space defined by this grid.
    @return cells.nfull() (for L0) or cells.nvertices() (for L1) */
    size_t ndata() const;

    void clear();

protected:
    void nc_read(netCDF::NcGroup *nc, std::string const &vname);
    void nc_write(netCDF::NcGroup *nc, std::string const &vname) const;
public:
    virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname);


    /** Remove cells and vertices not relevant to us --- for example, not in our MPI domain.
    This will be done AFTER we read it in.  It's an optimization. */
    void filter_cells(std::function<bool (long)> const &keep_fn);

};

class Grid_XY : public Grid
{
public:
    ~Grid_XY() {}

    /** Cell boundaries in the x direction.
    Sorted low to high.
    Number of grid cells in the x direction = x_boundaries.size() - 1. */
    std::vector<double> xb;

    /** Cell boundaries in the y direction.
    Sorted low to high.
    Number of grid cells in the y direction = y_boundaries.size() - 1. */
    std::vector<double> yb;

    int nx() const { return xb.size() - 1; }
    int ny() const { return yb.size() - 1; }

    void ncio(ibmisc::NcIO &ncio, std::string const &vname);
};

extern void sort_renumber_vertices(Grid &grid);
extern std::unique_ptr<Grid> new_grid(Grid::Type type);
extern std::unique_ptr<Grid> new_grid(ibmisc::NcIO &ncio, std::string const &vname);


class Grid_LonLat : public Grid
{
public:
    ~Grid_LonLat() {}

    /** Longitude of cell boundaries (degrees), sorted low to high.
    <b>NOTE:</b> lon_boundares.last() = 360.0 + lonb.first() */
    std::vector<double> lonb;

    /** Latitude of cell boundaries (degrees), sorted low to high.
    Does NOT include the polar "cap" cells.
    90 = north pole, -90 = south pole. */
    std::vector<double> latb;

    /** True if this grid contains a circular cap on the south pole.
    If not, then many triangular grid cells will meet at the south pole. */
    bool south_pole;

    /** True if this grid contains a circular cap on the north pole.
    If not, then many triangular grid cells will meet at the north pole. */
    bool north_pole;

    /** Number of segments used to represent each side of a grid cell
    with a polygonal approximation.  Increasing this number will
    decrease geometric error at the expense of computation time for
    the overlap matrix. */
    int points_in_side;

    /** Radius of Earth (m), or whatever planet you're looking at.
    Used to compute theoretical exact areas of graticules. That's different
    from the radius (or elliptical shape) used in the projection. */
    double eq_rad;

    /** Number of grid cell indices in longitude dimension */
    int nlon() const { return lonb.size() - 1; }

    /** Number of grid cell indices in latitude dimension */
    int nlat() const;

    void ncio(ibmisc::NcIO &ncio, std::string const &vname);
};  // class

}   // namespace

std::ostream &operator<<(std::ostream &out, icebin::Cell const &cell);
std::ostream &operator<<(std::ostream &out, icebin::Vertex const &vertex);

