#ifndef ICEBIN_GRIDSPEC_HPP
#define ICEBIN_GRIDSPEC_HPP

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
struct HntrSpec : public GridSpec {
    int im;    // Number of cells in east-west direction
    int jm;    // Number of cells in north-south direction

    // number (fraction) of cells in east-west direction from
    // International Date Line (180) to western edge of cell IA=1
    double offi;

    // minutes of latitude for non-polar cells on grid A
    double dlat;

    std::unique_ptr<GridSpec> clone()
        { return std::unique_ptr(new HntrSpec(*this)); }
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
    std::unique_ptr<HntrSpec> hntr;

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

}    // namespace
#endif
