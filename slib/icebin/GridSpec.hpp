#ifndef ICEBIN_GRIDSPEC_HPP
#define ICEBIN_GRIDSPEC_HPP

#include <boost/enum.hpp>
#include <ibmisc/netcdf.hpp>
#include <icebin/error.hpp>

namespace icebin {

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
    GridType type;

    //GridSpec() : type(GridType::GENERIC) {}
    GridSpec(GridType _type) : type(_type) {}
    virtual ~GridSpec() {}
    virtual long ncells_full() const = 0;
    virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname);
    virtual std::unique_ptr<GridSpec> clone() const = 0;
};
// -------------------------------------------------------------------
// ----------------------------------------------------
class GridSpec_Generic : public GridSpec {
    long _ncells_full;

public:
    GridSpec_Generic(long ncells_full = -1) : GridSpec(GridType::GENERIC), _ncells_full(ncells_full) {}

    long ncells_full() const { return _ncells_full; }
    void ncio(ibmisc::NcIO &ncio, std::string const &vname);
    std::unique_ptr<GridSpec> clone() const
        { return std::unique_ptr<GridSpec>(new GridSpec_Generic(*this)); }
};
// ----------------------------------------------------
struct GridSpec_XY : public GridSpec {

    /** Projection of this XY grid onto the sphere. */
    std::string sproj;

    /** Cell boundaries in the x direction.
    Sorted low to high.
    Number of grid cells in the x direction = x_boundaries.size() - 1. */
    std::vector<double> xb;

    /** Cell boundaries in the y direction.
    Sorted low to high.
    Number of grid cells in the y direction = y_boundaries.size() - 1. */
    std::vector<double> yb;

    /** (x,y) Dimensions in order of decreasing stride (same as Indexing class).
        (0,1) = (x,y) = x has largest stride
        (1,0) = (y,x) = y has largest stride */
    std::vector<int> indices;

    int nx() const { return xb.size() - 1; }
    int ny() const { return yb.size() - 1; }

    // -------------------------------------------
    long ncells_full() const { return nx() * ny(); }
    void ncio(ibmisc::NcIO &ncio, std::string const &vname);

    std::unique_ptr<GridSpec> clone() const
        { return std::unique_ptr<GridSpec>(new GridSpec_XY(*this)); }
    // -------------------------------------------


    // Used only when reading with ncio()
    GridSpec_XY() : GridSpec(GridType::XY) {}



    GridSpec_XY(
        std::string const &_sproj,
        std::vector<int> const &&_indices,
        std::vector<double> &&_xb,
        std::vector<double> &&_yb)
    : GridSpec(GridType::XY), sproj(_sproj), indices(std::move(_indices)),
        xb(std::move(_xb)),
        yb(std::move(_yb))
    {}



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
        std::string const &sproj,
        std::vector<int> const &&indices,    // Decreasing stride
        double x0, double x1, double dx,
        double y0, double y1, double dy);

    static GridSpec_XY make_with_centers(
        std::string const &sproj,
        std::vector<int> const &&indices,    // Decreasing stride
        double x0, double x1, double dx,
        double y0, double y1, double dy)
    {
        return make_with_boundaries(
            sproj, std::move(indices),
            x0-.5*dx, x1+.5*dx, dx,
            y0-.5*dy, y1+.5*dy, dy);
    }

};

// -------------------------------------------------------------------

/** Used to create a GridSpec_LonLat; not a GridSpec itself. */
struct HntrSpec {
    int im;    // Number of cells in east-west direction
    int jm;    // Number of cells in north-south direction

    // number (fraction) of cells in east-west direction from
    // International Date Line (180) to western edge of cell IA=1
    double offi;

    // minutes of latitude for non-polar cells on grid A
    double dlat;

    void ncio(ibmisc::NcIO &ncio, std::string const &vname);

    HntrSpec() : im(-1), jm(-1), offi(0.), dlat(0.) {}
    HntrSpec(int _im, int _jm, double _offi, double _dlat);

    int size() const { return im*jm; }
    int ndata() const { return size(); }    // Convention makes this more like regular ModelE grids

    bool is_set() const { return (im >= 0); }

    /** provides lon variable found in ModelE files, indicating
    longitude center of each cell. */
    std::vector<double> lonc() const;
    /** provides lat variable found in ModelE files, indicating
    latitude center of each cell. */
    std::vector<double> latc() const;
};


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

    /** (lon,lat) Dimensions in order of decreasing stride.
        (0,1) = lon has largest stride
        (1,0) = lat has largest stride (ModelE) */
    std::vector<int> indices;

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

    /** If this was created from a HntrSpec, here is that spec. */
    HntrSpec hntr;

    // --------------------------------------------
    GridSpec_LonLat() : GridSpec(GridType::LONLAT) {}
    GridSpec_LonLat(
        std::vector<double> &&_lonb,
        std::vector<double> &&_latb,
        std::vector<int> const &_indices,
        bool _south_pole,
        bool _north_pole,
        int _points_in_side,
        double _eq_rad,
        HntrSpec const &hntr = HntrSpec());

    void ncio(ibmisc::NcIO &ncio, std::string const &vname);

    // --------------------------------------------

    /** Number of grid cell indices in longitude dimension */
    int nlon() const { return lonb.size() - 1; }

    /** Number of grid cell indices in latitude dimension */
    int nlat() const;

    long ncells_full() const { return nlon() * nlat(); }


    /** @return [nlat()] latidue of cell centers */
    std::vector<double> latc() const;
    /** @return [nlon()] longitude of cell centers */
    std::vector<double> lonc() const;

    std::unique_ptr<GridSpec> clone() const
        { return std::unique_ptr<GridSpec>(new GridSpec_LonLat(*this)); }
};

/** Make a GridSpec_LonLat form a HntrSpec */
extern GridSpec_LonLat make_grid_spec(HntrSpec const &hntr, bool pole_caps, int points_in_side, double eq_rad);

extern void ncio_grid_spec(
    ibmisc::NcIO &ncio,
    std::unique_ptr<GridSpec> &spec,
    std::string const &vname);

}    // namespace
#endif
