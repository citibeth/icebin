#include <icebin/GridSpec.hpp>
#include <ibmisc/netcdf.hpp>

using namespace netCDF;
using namespace ibmisc;

namespace icebin {


void GridSpec::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    if (ncio.rw == 'w') {
        NcVar info_v = get_or_add_var(ncio, vname + ".info", "int", {});
        info_v.putAtt("type", std::string(type.str()));
    }
}


void GridSpec_Generic::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    GridSpec::ncio(ncio, vname);
    NcVar info_v = get_or_add_var(ncio, vname + ".info", "int", {});
    get_or_put_att(info_v, ncio.rw, "ncells_full", "int64", &_ncells_full, 1);
}



void GridSpec_XY::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    GridSpec::ncio(ncio, vname);

    auto xb_d = get_or_add_dim(ncio,
        vname + ".x_boundaries.length", this->xb.size());
    ncio_vector(ncio, this->xb, true,
        vname + ".x_boundaries", "double", {xb_d});

    auto yb_d = get_or_add_dim(ncio,
        vname + ".y_boundaries.length", this->yb.size());
    ncio_vector(ncio, this->yb, true,
        vname + ".y_boundaries", "double", {yb_d});

    NcVar info_v = get_or_add_var(ncio, vname + ".info", "int", {});
    get_or_put_att(info_v, ncio.rw, "sproj", sproj);
    if (ncio.rw == 'w') {
        int n;
        n = nx();
        get_or_put_att(info_v, ncio.rw, "nx", "int", &n, 1);
        n = ny();
        get_or_put_att(info_v, ncio.rw, "ny", "int", &n, 1);
    }
}
// ---------------------------------------------------------
GridSpec_LonLat::GridSpec_LonLat(
    std::vector<double> &&_lonb,
    std::vector<double> &&_latb,
    std::vector<int> const &_indices,
    bool _south_pole,
    bool _north_pole,
    int _points_in_side,
    double _eq_rad,
    HntrSpec const &_hntr)

: GridSpec(GridType::LONLAT), lonb(std::move(_lonb)), latb(std::move(_latb)),
indices(_indices),
south_pole(_south_pole), north_pole(_north_pole),
points_in_side(_points_in_side), eq_rad(_eq_rad), hntr(_hntr)
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
// ---------------------------------------------------------
GridSpec_LonLat make_grid_spec(HntrSpec const &hntr, bool pole_caps, int points_in_side, double eq_rad)
{
    if (hntr.im % 2 != 0) (*icebin_error)(-1,
        "IM must be even");

    std::vector<double> lonb;
    std::vector<double> latb;

    // Longitude grid boundaries
    double const deg_by_im = 360. / (double)hntr.im;
    for (int i=0; i<hntr.im; ++i) {
        // A longitude offset of -180 is added: in hntr grids,
        // the cell with longitude index=0 starts at 180 degrees W.
        // This splits the map in the Pacific Ocean when
        // plotting data naively
        lonb.push_back( -180. + (hntr.offi + (double)i) * deg_by_im );
    }
    lonb.push_back(lonb[0] + 360.);

    // Latitude grid boundaries
    latb.push_back(0);
    double dlat_d = hntr.dlat / 60.;    // Convert minutes -> degrees
    for (int j=1; j<hntr.jm/2; ++j) {
        double lat = j * dlat_d;
        latb.push_back(lat);
        latb.push_back(-lat);
    }
    if (!pole_caps) {
        double lat = hntr.jm/2 * dlat_d;
        if (std::abs(lat-90.) < 1.e-10) lat = 90.;
        latb.push_back(lat);
        latb.push_back(-lat);
    }

    std::sort(latb.begin(), latb.end());

    GridSpec_LonLat spec(
        std::move(lonb), std::move(latb),
        {1,0},    // Latitude has largest stride
        pole_caps, pole_caps,
        points_in_side,
        eq_rad, hntr);
    return spec;
}
// -----------------------------------------------------
HntrSpec::HntrSpec(int _im, int _jm, double _offi, double _dlat)
    : im(_im), jm(_jm), offi(_offi), dlat(_dlat)
{
    // Check for common error of degrees instead of minutes
    // (this is heuristic)
    if (dlat < .1 * (180.*60./jm)) (*icebin_error)(-1,
        "dlat in HntrGrid(%d,%d,%f,%f) seems to small; it should be in MINUTES ,not DEGREES: %f", im,jm,offi,dlat, .1 * (180.*60./jm));
}
// -----------------------------------------------------
void HntrSpec::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    auto hntr_v = get_or_add_var(ncio, vname, "int64", {});
    if (ncio.rw == 'w') hntr_v.putAtt("comment",
        "Parameters to instantiate a icebin::modele::HntrGrid description of this grid");

    get_or_put_att(hntr_v, ncio.rw, "im", "int", &im, 1);
    get_or_put_att(hntr_v, ncio.rw, "jm", "int", &jm, 1);
    get_or_put_att(hntr_v, ncio.rw, "offi", "double", &offi, 1);
    get_or_put_att(hntr_v, ncio.rw, "dlat", "double", &dlat, 1);
}

void GridSpec_LonLat::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    GridSpec::ncio(ncio, vname);

    // Read the HntrGrid, if that's how we were made
    std::string const hntr_vname = vname + ".hntr";
    if (ncio.rw == 'r') {
        auto hntr_v = ncio.nc->getVar(hntr_vname);
        if (!hntr_v.isNull()) {
            hntr.ncio(ncio, hntr_vname);
        }
    } else {
        if (hntr.is_set()) {
            hntr.ncio(ncio, hntr_vname);
        }
    }

    auto lonb_d = get_or_add_dim(ncio,
        vname + ".lon_boundaries.length", this->lonb.size());
    ncio_vector(ncio, this->lonb, true,
        vname + ".lon_boundaries", "double", {lonb_d});

    auto latb_d = get_or_add_dim(ncio,
        vname + ".lat_boundaries.length", this->latb.size());
    ncio_vector(ncio, this->latb, true,
        vname + ".lat_boundaries", "double", {latb_d});

    NcVar info_v = get_or_add_var(ncio, vname + ".info", "int", {});

    get_or_put_att(info_v, ncio.rw, "eq_rad", "double", &eq_rad, 1);
    if (ncio.rw == 'w') info_v.putAtt("eq_rad.comment",
        "Radius of Earth (m), or whatever planet you're looking at.  "
        "Used to compute theoretical exact areas of graticules. That's different "
        "from the radius (or elliptical shape) used in the projection.");

    get_or_put_att(info_v, ncio.rw, "north_pole_cap", north_pole);
    get_or_put_att(info_v, ncio.rw, "south_pole_cap", south_pole);
    get_or_put_att(info_v, ncio.rw, "points_in_side", "int", &points_in_side, 1);
    if (ncio.rw == 'w') {
        int n;
        n = nlon();
        get_or_put_att(info_v, ncio.rw, "nlon", "int", &n, 1);
        n = nlat();
        get_or_put_att(info_v, ncio.rw, "nlat", "int", &n, 1);
    }
}
// ---------------------------------------------------------
/** provides lon variable found in ModelE files, indicating
longitude center of each cell. */
std::vector<double> HntrSpec::lonc() const
{
    std::vector<double> ret;

    // Longitude grid boundaries
    double const deg_by_im = 360. / (double)im;
    for (int i=0; i<im; ++i) {
        // A longitude offset of -180 is added: in hntr grids,
        // the cell with longitude index=0 starts at 180 degrees W.
        // This splits the map in the Pacific Ocean when
        // plotting data naively
        ret.push_back( -180. + (offi + (double)i + .5) * deg_by_im );
    }
    return ret;
}

/** provides lat variable found in ModelE files, indicating
latitude center of each cell. */
std::vector<double> HntrSpec::latc() const
{
    std::vector<double> ret;

    // Latitude grid boundaries
    double dlat_d = dlat / 60.;    // Convert minutes -> degrees
    for (int j=0; j<jm/2; ++j) {
        double lat = ((double)j + .5) * dlat_d;
        ret.push_back(lat);
        ret.push_back(-lat);
    }

    std::sort(ret.begin(), ret.end());
    return ret;
}

// ---------------------------------------------------------
int GridSpec_LonLat::nlat() const {
    int south_pole_offset, north_pole_offset;

    // Get around GCC (4.9.3) bug when converting uninitialized bools
    south_pole_offset = (south_pole ? 2 : 0);
    north_pole_offset = (north_pole ? 2 : 0);

    south_pole_offset >>= 1;
    north_pole_offset >>= 1;

    return latb.size() - 1 + south_pole_offset + north_pole_offset;
}

// -------------------------------------------------
void ncio_grid_spec(
    NcIO &ncio,
    std::unique_ptr<GridSpec> &spec,
    std::string const &vname)
{
    if (ncio.rw == 'w') {
        spec->ncio(ncio, vname);
    } else {
        NcVar info_v = get_or_add_var(ncio, vname + ".info", "int", {});
        GridType type;
        get_or_put_att_enum(info_v, ncio.rw, "type", type);

        switch(type.index()) {
            case GridType::GENERIC :
                spec.reset(new GridSpec_Generic());
            break;
            case GridType::XY :
                spec.reset(new GridSpec_XY());
            break;
            case GridType::LONLAT :
                spec.reset(new GridSpec_LonLat());
            break;
            default :
                (*icebin_error)(-1,
                    "Unrecognized GridType: %s", type.str());
        }
        spec->ncio(ncio, vname);
    }
}

}    // namespace

