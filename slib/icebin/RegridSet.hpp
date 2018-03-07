#ifndef ICEBIN_REGRID_SET_HPP
#define ICEBIN_REGRID_SET_HPP

#include <unordered_set>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>

#include <icebin/eigen_types.hpp>

namespace icebin {

class IceRegridder;


// -----------------------------------------------------------
class RegridSet {

    /** All matrices in the RegridSet share these params */
    RegridParams _params;

public:
    RegridParams const &params() { return _params; }

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
    virtual std::unique_ptr<ibmisc::lintransform::Weighted> matrix_s(
        std::string const &spec_name) const;
};
// -----------------------------------------------------------------
}    // namespace
#endif    // guard
