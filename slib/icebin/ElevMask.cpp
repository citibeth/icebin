#include <icebin/ElevMask.hpp>
#include <ibmisc/blitz.hpp>

using namespace ibmisc;

namespace icebin {

// Types of PISM cover
// From pism/src/base/util/Mask.hh
struct IceMask {
    static const char UNKNOWN          = -1;
    static const char ICE_FREE_BEDROCK = 0;
    static const char GROUNDED_ICE     = 2;
    static const char FLOATING_ICE     = 3;
    static const char ICE_FREE_OCEAN   = 4;
};

static double const NaN = std::numeric_limits<double>::quiet_NaN();

/** Reads and allocate elevmaskI arrays from a PISM state file.
@param emI Elevation-mask for continental area (ice and bare land)
@param emI_ice Elevation-mask for just ice-covered areas */
void read_elevmask_pism(
    std::string const &fname,
    blitz::Array<double,1> &emI_land,
    blitz::Array<double,1> &emI_ice)
{
    // Variables to be used from the PISM file.
    blitz::Array<double,2> topg, thk;
    blitz::Array<char,2> mask;

    /* Read it in */
    {NcIO ncio(fname, 'r');
        ncio_blitz_alloc(ncio, topg, "topg", "double");
        ncio_blitz_alloc(ncio, thk, "thk", "double");
        ncio_blitz_alloc(ncio, mask, "mask", "");
    }

    // Move to 1D
    auto topg1(reshape1(topg));
    auto thk1(reshape1(thk));
    auto mask1(reshape1(mask));
    auto nI(topg1.extent(0));

    /* Generate emI_land and emI_ice */
    emI_land.reference(blitz::Array<double,1>(nI));
    emI_ice.reference(blitz::Array<double,1>(nI));

    for (int iI=0; iI<nI; ++iI) {
        switch(mask1(iI)) {
            case IceMask::GROUNDED_ICE :
            case IceMask::FLOATING_ICE :
                emI_land(iI) = emI_ice(iI) = topg(iI) + thk(iI);
            break;
            case IceMask::ICE_FREE_OCEAN :
            case IceMask::UNKNOWN :
                emI_land(iI) = emI_ice(iI) = NaN;
            break;
            case IceMask::ICE_FREE_BEDROCK :
                emI_land(iI) = topg(iI);
                emI_ice(iI) = NaN;
            break;
        }
    }
}

}    // namespace
