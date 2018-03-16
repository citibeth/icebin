#ifndef ICEBIN_REGRID_MATRICES_HPP
#define ICEBIN_REGRID_MATRICES_HPP

#include <unordered_set>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <ibmisc/linear/linear.hpp>

#include <icebin/eigen_types.hpp>

namespace icebin {

class IceRegridder;


/** Parameters controlling the generation of regridding matrices */
struct RegridParams {
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

    RegridParams() : scale(true), correctA(false), sigma({0.,0.,0.}) {}

    RegridParams(bool _scale, bool _correctA, std::array<double,3> const &_sigma) :
        scale(_scale), correctA(_correctA), sigma(_sigma) {}
};
// -----------------------------------------------------------
class RegridMatrices {

    /** All matrices in the RegridMatrices share these params */
    RegridParams _params;

public:
    RegridMatrices(RegridParams const &params) : _params(params) {}
    virtual ~RegridMatrices() {}

    RegridParams const &params() const { return _params; }

    /** Retrieves a regrid matrix, and weights (area) of the input and
        output grid cells.
    @param spec_name: The matrix to produce.
        Should be "AvI", "IvA", "EvI", "IvE", "AvE" or "EvA".
    @param scale: Produce scaled matrix?
        true  --> [kg m-2]
        false --> [kg]
    @param correctA: Correct for projection error in A or E grids?
    @return The regrid matrix and weights
    */
    virtual std::unique_ptr<ibmisc::linear::Weighted> matrix(
        std::string const &spec_name) const = 0;
};
// -----------------------------------------------------------------
}    // namespace
#endif    // guard
