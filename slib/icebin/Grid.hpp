/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
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
#include <icebin/GridSpec.hpp>

#ifdef BUILD_MODELE
#include <icebin/modele/hntr.hpp>
#endif

namespace icebin {

class Cell;

// --------------------------------------------------
struct Point {
    double const x;
    double const y;

    Point() : x(0), y(0) {}
    Point(double _x, double _y) : x(_x), y(_y) {}
};

struct Vertex : public Point {
    long index;

    Vertex() : index(-1) {}

    Vertex(double _x, double _y, int _index=-1) :
        Point(_x,_y), index(_index) {}

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

    double proj_area(ibmisc::Proj_LL2XY const *proj) const;   // OPTIONAL

    Point centroid() const;
};      // class Cell
// ----------------------------------------------------
class Grid;

class GridGen {};  // Tagging class

/** Specialized dict-like structure used for cells and vertices in a grid. */
template<class CellT>
class GridMap {
    friend class Grid;
protected:
    typedef std::unordered_map<long, std::unique_ptr<CellT>> MapT;
    MapT _cells;
    long _nfull = -1;
    long _max_realized_index = -1;

public:
    GridMap(long nfull) : _nfull(nfull) {}
    GridMap() : _nfull(-1) {}

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



// -------------------------------------------------------------------
class Grid {
    // These would all be const, except for ncio()
    // Instead, we use std::unique_ptr<const Grid> to const-ify everything
public:
    GridMap<Vertex> vertices;
    GridMap<Cell> cells;

    // GridType type;    // now spec->type
    std::unique_ptr<GridSpec> spec;    // Details used to generate the grid
    GridCoordinates coordinates;
    GridParameterization parameterization;

    /** Conversion between n-dimensional indexing used natively on the
    grid, and 1-D indexing used in IceBin. */
    ibmisc::Indexing indexing;

    std::string name;

    /** If scoord == "xy": The projection that relates x,y coordinates
    here to a specific point on the globe (as a Proj.4 String). */
    std::string sproj;

    // Just used for ncio() read
    Grid() {}

    Grid(
        std::string const &_name,
        std::unique_ptr<GridSpec> &&_spec,
        GridCoordinates _coordinates,
        std::string const &_sproj,
        GridParameterization _parameterization,

        /** Conversion between n-dimensional indexing used natively on the
        grid, and 1-D indexing used in IceBin. */
        ibmisc::Indexing &&_indexing,

        GridMap<Vertex> &&_vertices,
        GridMap<Cell> &&_cells);


    /** The size of the vector space defined by this grid.
    @return cells.nfull() (for L0) or cells.nvertices() (for L1) */
    size_t ndata() const;

    /** The number of grid cells / basis vectors used in the vector
    space. */
    size_t nrealized() const;

    void clear();

    /** For now, just return the geographic center of the cell's polygon.
        But this might be revisited for finite element */
    Point centroid(Cell const &cell) const
        { return cell.centroid(); }

protected:
    void nc_read(netCDF::NcGroup *nc, std::string const &vname);
    void nc_write(netCDF::NcGroup *nc, std::string const &vname) const;
public:
    virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname, bool rw_full=true);


    /** Remove cells and vertices not relevant to us --- for example, not in our MPI domain.
    This will be done AFTER we read it in.  It's an optimization. */
    void filter_cells(std::function<bool (long)> const &keep_fn);
};

void sort_renumber_vertices(Grid &grid);


}   // namespace

std::ostream &operator<<(std::ostream &out, icebin::Cell const &cell);
std::ostream &operator<<(std::ostream &out, icebin::Vertex const &vertex);

