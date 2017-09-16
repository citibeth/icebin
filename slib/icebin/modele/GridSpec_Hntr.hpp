#ifndef GRIDSPECHNTR_HPP
#define GRIDSPECHNTR_HPP

#include <ibmisc/indexing.hpp>
#include <icebin/Grid.hpp>
#include <icebin/modele/hntr.hpp>

namespace icebin {
namespace modele {

struct GridSpec_Hntr {
    std::string name;
    HntrGrid hntr;
    std::function<bool(long, double, double, double, double)> spherical_clip;

    /** True if this grid contains a circular cap on the north pole.
    If not, then many triangular grid cells will meet at the north pole. */
    bool pole_caps = false;

    /** Number of segments used to represent each side of a grid cell
    with a polygonal approximation.  Increasing this number will
    decrease geometric error at the expense of computation time for
    the overlap matrix. */
    int points_in_side = 1;

    /** Radius of Earth (m), or whatever planet you're looking at.
    Used to compute theoretical exact areas of graticules. That's different
    from the radius (or elliptical shape) used in the projection. */
    double eq_rad = 1.;


    GridSpec_Hntr(HntrGrid const &_hntr);

    void make_grid(Grid_LonLat &grid);
};



}}    // namespace

#endif    // guard
