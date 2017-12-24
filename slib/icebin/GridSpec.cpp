#include <icebin/GridSpec.hpp>

namespace icebin {

void GridSpec_XY::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{

    auto xb_d = get_or_add_dim(ncio,
        vname + ".x_boundaries.length", this->xb.size());
    ncio_vector(ncio, this->xb, true,
        vname + ".x_boundaries", "double", {xb_d});

    auto yb_d = get_or_add_dim(ncio,
        vname + ".y_boundaries.length", this->yb.size());
    ncio_vector(ncio, this->yb, true,
        vname + ".y_boundaries", "double", {yb_d});

    NcVar info_v = get_or_add_var(ncio, vname + ".info", "int", {});
    if (ncio.rw == 'w') {
        int n;
        n = nx();
        get_or_put_att(info_v, ncio.rw, "nx", "int", &n, 1);
        n = ny();
        get_or_put_att(info_v, ncio.rw, "ny", "int", &n, 1);
    }
}
// ---------------------------------------------------------


GridSpec_Lonlat::GridSpec_LonLat(GridSpec_Hntr const &_hntr)
{
    hntr.reset(new GridSpec_Hntr(hntr));

    if (hntr.im % 2 != 0) (*icebin_error)(-1,
        "IM must be even");

    GridGen_LonLat spec;
    spec.name = name;
    spec.spherical_clip = spherical_clip;
    spec.points_in_side = points_in_side;
    spec.eq_rad = eq_rad;

    // Longitude grid boundaries
    double const deg_by_im = 360. / (double)hntr.im;
    for (int i=0; i<hntr.im; ++i) {
        spec.lonb.push_back( 180. + (hntr.offi + (double)i) * deg_by_im );
    }
    spec.lonb.push_back(spec.lonb[0] + 360.);

    // Latitude grid boundaries
    spec.latb.push_back(0);
    double dlat_d = hntr.dlat / 60.;    // Convert minutes -> degrees
    for (int j=1; j<hntr.jm/2; ++j) {
        double lat = j * dlat_d;
        spec.latb.push_back(lat);
        spec.latb.push_back(-lat);
    }
    if (pole_caps) {
        spec.north_pole = true;
        spec.south_pole = true;

    } else {
        spec.north_pole = false;
        spec.south_pole = false;

        double lat = hntr.jm/2 * dlat_d;
        if (std::abs(lat-90.) < 1.e-10) lat = 90.;
        spec.latb.push_back(lat);
        spec.latb.push_back(-lat);
    }

    std::sort(spec.latb.begin(), spec.latb.end());

    spec.indexing = Indexing(
        {"lon", "lat"}, {0,0}, {spec.nlon(), spec.nlat()}, {1,0});  // col major
    spec.make_grid(grid);




}


void GridSpec_LonLat::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
#ifdef BUILD_MODELE
    // Read the HntrGrid, if that's how we were made
    if (ncio.rw == 'r') {
        std::string const hntr_vname = vname + ".hntr";
        auto hntr_v = ncio.nc->getVar(hntr_vname);
        if (!hntr_v.isNull()) {
            HntrGrid _hntr;
            _hntr.ncio(ncio, hntr_vname);
            hntr.reset(new HntrGrid(std::move(_hntr)));
        }
    } else {
        if (hntr.get()) {
            HntrGrid _hntr(*hntr);
            _hntr.ncio(ncio, vname + ".hntr");
        }
    }
#endif

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

int GridSpec_LonLat::nlat() const {
    int south_pole_offset, north_pole_offset;

    // Get around GCC bug when converting uninitialized bools
    south_pole_offset = (south_pole ? 2 : 0);
    north_pole_offset = (north_pole ? 2 : 0);

    south_pole_offset >>= 1;
    north_pole_offset >>= 1;

    return latb.size() - 1 + south_pole_offset + north_pole_offset;
}

}    // namespace
