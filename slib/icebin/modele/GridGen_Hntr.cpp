#include <icebin/modele/GridGen_Hntr.hpp>
#include <icebin/Grid.hpp>
#include <ibmisc/indexing.hpp>
#include <icebin/gridgen/GridGen_LonLat.hpp>

using namespace ibmisc;
using namespace std::placeholders;  // for _1, _2, _3...

namespace icebin {
namespace modele {

GridGen_Hntr::GridGen_Hntr(HntrGrid const &_hntr) : hntr(_hntr) {}

void GridGen_Hntr::make_grid(Grid_LonLat &grid)
{
    grid.hntr.reset(new HntrGrid(hntr));

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



}}
