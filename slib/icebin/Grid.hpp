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

    double proj_area(ibmisc::Proj_LL2XY const *proj);   // OPTIONAL

    Point centroid() const;
};      // class Cell
// ----------------------------------------------------
class Grid;

class GridGen {};  // Tagging class

/** Specialized dict-like structure used for cells and vertices in a grid. */
template<class CellT>
class GridMap {
//    friend class Grid;
//    friend class GridGen_XY;
//    friend class GridGen_LonLat;
//    friend class GridGen_Exchange;
protected:
    typedef std::unordered_map<long, std::unique_ptr<CellT>> MapT;
    MapT _cells;
    long _nfull = -1;
    long _max_realized_index = -1;

    GridMap(long nfull) : _nfull(nfull) {}

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
// Corresponds to classes
BOOST_ENUM_VALUES( GridType, int,
    (GENERIC)       (0)     // Just use the Grid base class
    (XY)            (1)     // Rectilinear X/Y grid
    (LONLAT)        (2)     // Global lat-lon grid (maybe with polar caps)
//    (EXCHANGE)      (3)     // Exchange grid, from overlap of two other grids
    (CUBESPHERE)    (4)     // Global Cubed Sphere grid
    (MESH)          (5)     // Arbitrary mesh (could be global or on plane)
)

BOOST_ENUM_VALUES( GridCoordinates, int,
    (XY)            (0)     // Vertices in x/y coordinates on a plane
    (LONLAT)        (1)     // Vertices in lon/lat coordinates on a sphere
)

BOOST_ENUM_VALUES( GridParameterization, int,
    (L0)            (0)     // Constant value in each grid cell
    (L1)            (1)     // Value specified at each vertex, slope inbetween
)
// -------------------------------------------------------------------
class GridSpec
{
public:
    const GridType type;

    GridSpec() : type(GridType::GENERIC) {}
    GridSpec(GridType _type) : type(_type) {}
    virtual ~GridSpec() {}
    virtual long ncells_full() = 0;
    virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname) = 0;
    virtual std::unique_ptr<GridSpec> clone() = 0;
};
// ----------------------------------------------------
struct GridSpec_XY : public GridSpec {

    /** Cell boundaries in the x direction.
    Sorted low to high.
    Number of grid cells in the x direction = x_boundaries.size() - 1. */
    std::vector<double> const xb;

    /** Cell boundaries in the y direction.
    Sorted low to high.
    Number of grid cells in the y direction = y_boundaries.size() - 1. */
    std::vector<double> const yb;

    int nx() const { return xb.size() - 1; }
    int ny() const { return yb.size() - 1; }

    GridSpec_XY(
        std::vector<double> &&_xb,
        std::vector<double> &&_yb)
    : GridSpec(GridType::XY), xb(std::move(_xb)), yb(std::move(_yb)) {}


    /** Create a new Cartesian grid with evenly spaced grid cell boundaries.
    @param name Value of <generic-name>.info:name in netCDF file.
    @param x0 Lowest boundary in the x direction.
    @param x1 Highest boundary in the x direction.  
    @param dx Size of grid cell in the x direction.
        Will be adjusted if (x1-x0) is not an even multiple of dx
    @param y0 Lowest boundary in the y direction.
    @param y1 Highest boundary in the y direction.  
    @param dy Size of grid cell in the y direction.
        Will be adjusted if (y1-y0) is not an even multiple of dy
    @param euclidian_clip Only realize grid cells that pass this test.
    @return The newly created GridGen_XY.
    @see EuclidianClip
    */
    static GridSpec_XY make_with_boundaries(
        double x0, double x1, double dx,
        double y0, double y1, double dy);

    static GridSpec_XY make_with_centers(
        double x0, double x1, double dx,
        double y0, double y1, double dy)
    {
        return GridSpec_XY::set_xy_boundaries(spec,
            x0-.5*dx, x1+.5*dx, dx,
            y0-.5*dy, y1+.5*dy, dy);
    }

    std::unique_ptr<GridSpec> clone()
        { return std::unique_ptr(new GridSpec_XY(*this)); }
};
// -------------------------------------------------------------------
/** Represents a lat-lon grid on a sphere with non-equally-spaced grid
cell boundaries.  Can optionally represent circular "cap" cells at the
north and south pole, if needed. */
struct GridSpec_LonLat : public GridSpec {

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

    /** Number of segments (not points) used to represent each side of a grid cell
    with a polygonal approximation.  Increasing this number will
    decrease geometric error at the expense of computation time for
    the overlap matrix.  Must be at least 1. */
    int points_in_side;

    /** Radius of Earth (m), or whatever planet you're looking at.
    Used to compute theoretical exact areas of graticules. That's different
    from the radius (or elliptical shape) used in the projection. */
    double eq_rad;

    // --------------------------------------------
    GridSpec_LonLat(
        std::vector<double> &&_lonb,
        std::vector<double> &&_latb,
        bool _south_pole,
        bool _north_pole,
        int _points_in_side,
        double _eq_rad)
    : GridSpec(GridType::LONLAT), lonb(std::move(_lonb)), latb(std::move(_latb)),
    south_pole(_south_pole), north_pole(_north_pole),
    points_in_side(_points_in_side), eq_rad(_eq_rad)
    {
        // Error-check the input parameters
        if (south_pole && latb[0] == -90.0) {
            (*icebin_error)(-1,
                "latb[] cannot include -90.0 if you're including the south pole cap");
        }
        if (north_pole && latb.back() == 90.0) {
            (*icebin_error)(-1,
                "latb[] cannot include 90.0 if you're including the north pole cap");
        }
    }

    // --------------------------------------------

    /** Number of grid cell indices in longitude dimension */
    int nlon() const { return lonb.size() - 1; }

    /** Number of grid cell indices in latitude dimension */
    int nlat() const {
        const int south_pole_offset = (south_pole ? 1 : 0);
        const int north_pole_offset = (north_pole ? 1 : 0);
        return latb.size() - 1 + south_pole_offset + north_pole_offset;
    }

    /** @return [nlat()] latidue of cell centers */
    std::vector<double> latc() const;
    /** @return [nlon()] longitude of cell centers */
    std::vector<double> lonc() const;

    std::unique_ptr<GridSpec> clone()
        { return std::unique_ptr(new GridSpec_LonLat(*this)); }
};

struct GridSpec_Hntr : public GridSpec {
    int im;    // Number of cells in east-west direction
    int jm;    // Number of cells in north-south direction

    // number (fraction) of cells in east-west direction from
    // International Date Line (180) to western edge of cell IA=1
    double offi;

    // minutes of latitude for non-polar cells on grid A
    double dlat;

    std::unique_ptr<GridSpec> clone()
        { return std::unique_ptr(new GridSpec_Hntr(*this)); }
};

// -------------------------------------------------------------------
class Grid {
    // These would all be const, except for ncio()
    // Instead, we use std::unique_ptr<const Grid> to const-ify everything
public:
    GridMap<Vertex> vertices;
    GridMap<Cell> cells;

    Type type;
    std::unique_ptr<GridSpec> spec;    // Details used to generate the grid
    Coordinates coordinates;
    Parameterization parameterization;

    /** Conversion between n-dimensional indexing used natively on the
    grid, and 1-D indexing used in IceBin. */
    ibmisc::Indexing indexing;

    std::string name;

    /** If scoord == "xy": The projection that relates x,y coordinates
    here to a specific point on the globe (as a Proj.4 String). */
    std::string sproj;

    // Just used for ncio() read
    Grid(GridType _type)
    : type(_type)
    {
        switch(type.index()) {
            case Grid::Type::GENERIC :
                spec = std::unique_ptr<GridSpec>(new GridSpec(GridType::GENERIC));
            case Grid::Type::XY :
                spec = std::unique_ptr<GridSpec>(new GridSpec_XY());
            case Grid::Type::LONLAT :
                spec = std::unique_ptr<GridSpec>(new GridSpec_LonLat());
            default :
                (*icebin_error)(-1,
                    "Unrecognized Grid::Type: %s", type.str());
        }
    }

    Grid(
        std::string const &_name,
        GridType _type,
        GridCoordinates _coordinates,
        std::string const &_sproj,
        GridParameterization _parameterization,

        /** Conversion between n-dimensional indexing used natively on the
        grid, and 1-D indexing used in IceBin. */
        ibmisc::Indexing const &_indexing,

        std::unique_ptr<GridSpec> &&_spec,
        GridMap<Vertex> &&_vertices,
        GridMap<Cell> &&_cells)
    : name(_name), type(_type), spec(std::move(_spec)),
    coordinates(_coordinates), sproj(_sproj), parameterization(_parameterization),
    indexing(_indexing), vertices(std::move(_vertices)), cells(std::move(_cells))
    {}



    /** The size of the vector space defined by this grid.
    @return cells.nfull() (for L0) or cells.nvertices() (for L1) */
    size_t ndata() const;

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





class Grid_XY : public GridExtra
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

    void ncio(ibmisc::NcIO &ncio, std::string const &vname, bool rw_full=true);
};

extern void sort_renumber_vertices(Grid &grid);
extern std::unique_ptr<Grid> new_grid(ibmisc::NcIO &ncio, std::string const &vname);


class Grid_LonLat : public GridExtra
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

#ifdef BUILD_MODELE
    /** The grid description, as needed for Hntr regridding algorithm. */
    std::unique_ptr<modele::HntrGrid const> hntr;
#endif


    /** Number of grid cell indices in longitude dimension */
    int nlon() const { return lonb.size() - 1; }

    /** Number of grid cell indices in latitude dimension */
    int nlat() const;

    void ncio(ibmisc::NcIO &ncio, std::string const &vname, bool rw_full=true);
};  // class Grid_LonLat

}   // namespace

std::ostream &operator<<(std::ostream &out, icebin::Cell const &cell);
std::ostream &operator<<(std::ostream &out, icebin::Vertex const &vertex);

