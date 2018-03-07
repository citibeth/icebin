#ifndef ICEBIN_REGRID_MATRICES_DYNAMIC_HPP
#define ICEBIN_REGRID_MATRICES_DYNAMIC_HPP

#include <unordered_set>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>

#include <icebin/eigen_types.hpp>

namespace icebin {

class IceRegridder;


// -----------------------------------------------------------
/** Holds the set of "Ur" (original) matrices produced by an
    IceRegridder for a SINGLE ice sheet. */
class RegridMatrices_Dynamic : public RegridMatrices {
public:
    /** ice_regridder from which this was constructed */
    IceRegridder const * const ice_regridder;

    ibmisc::TmpAlloc tmp;    // Stores local vars for different types.  TODO: Maybe re-do this as simple classmember variables.  At least, see where it used (by removing it and running the compiler)

    typedef std::function<std::unique_ptr<ibmisc::lintransform::Weighted>(
        std::array<SparseSetT *,2> dims, RegridMatrices::Params const &params)> MatrixFunction;

    std::map<std::string, MatrixFunction> regrids;

    RegridMatrices(IceRegridder const * const _ice_regridder) : ice_regridder(_ice_regridder) {}

    void add_regrid(std::string const &spec,
        MatrixFunction const &regrid);


    /** Retrieves a regrid matrix, and weights (area) of the input and
        output grid cells.
    @param spec_name: The matrix to produce.
        Should be "AvI", "IvA", "EvI", "IvE", "AvE" or "EvA".
    @param scale: Produce scaled matrix?
        true  --> [kg m-2]
        false --> [kg]
    @param correctA: Correct for projection error in A or E grids?
    @return The regrid matrix and weights.  All in dense indexing.
    */
    std::unique_ptr<ibmisc::lintransform::Eigen> matrix_d(
        std::string const &spec_name,
        std::array<SparseSetT *,2> dims,
        RegridMatrices::Params const &_params) const;

    // Virtual function promises a matrix in sparse indexing

};
// -----------------------------------------------------------------
}    // namespace
#endif    // guard
