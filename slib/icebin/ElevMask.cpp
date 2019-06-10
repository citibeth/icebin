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
    int const itime,    // Value of PISM time dimension to read
    blitz::Array<double,1> &emI_land,
    blitz::Array<double,1> &emI_ice)
{
    {NcIO ncio(fname, 'r');
        netCDF::NcVar nc_topg = ncio.nc->getVar("topg");
        std::vector<NamedDim> ndims(named_dims(nc_topg));
    }

    // Read the file
    blitz::Array<double,2> topg(ndims[1].extent, ndims[2].extent);
    blitz::Array<double,2> thk(ndims[1].extent, ndims[2].extent);
    blitz::Array<int8_t,2> mask(ndims[1].extent, ndims[2].extent);

    /* Read it in */
    {NcIO ncio(fname, 'r');
        ncio_blitz_partial(ncio, topg, "topg", "double", {}, {itime,0,0}, {1,2});
        ncio_blitz_partial(ncio, thk, "thk", "double", {}, {itime,0,0}, {1,2});
        ncio_blitz_partial(ncio, mask, "mask", "double", {}, {itime,0,0}, {1,2});
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
                emI_land(iI) = emI_ice(iI) = topg1(iI) + thk1(iI);
            break;
            case IceMask::ICE_FREE_OCEAN :
            case IceMask::UNKNOWN :
                emI_land(iI) = emI_ice(iI) = NaN;
            break;
            case IceMask::ICE_FREE_BEDROCK :
                emI_land(iI) = topg1(iI);
                emI_ice(iI) = NaN;
            break;
        }
    }
}

}    // namespace
