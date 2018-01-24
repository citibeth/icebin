#include <icebin/modele/GridGen_Hntr.hpp>
#include <icebin/Grid.hpp>
#include <ibmisc/indexing.hpp>
#include <icebin/gridgen/GridGen_LonLat.hpp>

using namespace ibmisc;
using namespace std::placeholders;  // for _1, _2, _3...

namespace icebin {
namespace modele {

void GridSpec_LonLat make_grid_spec(HntrSpec &hntr, bool pole_caps, int points_in_side, double eq_rad)
{
    if (hntr.im % 2 != 0) (*icebin_error)(-1,
        "IM must be even");

    std::vector<double> lonb;
    std::vector<double> latb;

    // Longitude grid boundaries
    double const deg_by_im = 360. / (double)hntr.im;
    for (int i=0; i<hntr.im; ++i) {
        lonb.push_back( 180. + (hntr.offi + (double)i) * deg_by_im );
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
        eq_rad);

    spec.hntr.reset(new HntrSpec(hntr));
    return spec;
}



}}
