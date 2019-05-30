#ifndef ICEBIN_ELEVMASK_HPP
#define ICEBIN_ELEVMASK_HPP

#include <blitz/array.h>

namespace icebin {

// From pism/src/base/util/Mask.hh
struct IceMask {
    static const char UNKNOWN          = -1;
    static const char ICE_FREE_BEDROCK = 0;
    static const char GROUNDED_ICE     = 2;
    static const char FLOATING_ICE     = 3;
    static const char ICE_FREE_OCEAN   = 4;
};

template<int RANK>
struct ElevMask {
    blitz::Array<double,RANK> elev;
    blitz::Array<char,RANK> mask;

    ElevMask(
        blitz::Array<double,RANK> const &_elev,
        blitz::Array<char,RANK> const &_mask)
    : elev(_elev), mask(_mask) {}
};

}    // namespace
#endif    // guard
