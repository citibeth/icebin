#ifndef ICEBIN_REGRIDPARAMS_HPP
#define ICEBIN_REGRIDPARAMS_HPP

#include <array>

namespace icebin {

    /** Parameters controlling the generation of regridding matrices */
    struct RegridMatrixParams {
        /** Produce a scaled vs. unscaled matrix (Scaled means divide
        by the weights vector).  Used for all matrices. */
        bool scale;

        /** Correct for changes in area due to projections?  Used for
        all matrices. */
        bool correctA;

        /** If non-zero, smooth the resulting matrix, with sigma as the
        scale length of the smoothing.  Used for IvA and IvE. */
        std::array<double,3> sigma;

        /** Tells if these parameters are asking us to smooth */
        bool smooth() const { return sigma[0] != 0; }

        Params() : scale(true), correctA(false), sigma({0.,0.,0.}) {}

        Params(bool _scale, bool _correctA, std::array<double,3> const &_sigma) :
            scale(_scale), correctA(_correctA), sigma(_sigma) {}
    };

}    // namespace
#endif
