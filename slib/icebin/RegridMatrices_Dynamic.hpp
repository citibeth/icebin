#ifndef ICEBIN_REGRID_MATRICES_DYNAMIC_HPP
#define ICEBIN_REGRID_MATRICES_DYNAMIC_HPP

#include <unordered_set>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <ibmisc/linear/eigen.hpp>
#include <icebin/RegridMatrices.hpp>
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

    typedef std::function<std::unique_ptr<ibmisc::linear::Weighted_Eigen>(
        std::array<SparseSetT *,2> dims, RegridParams const &params)> MatrixFunction;

    std::map<std::string, MatrixFunction> regrids;

    RegridMatrices_Dynamic(
        IceRegridder const * const _ice_regridder,
        RegridParams const &params)
    : RegridMatrices(params), ice_regridder(_ice_regridder) {}

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
    std::unique_ptr<ibmisc::linear::Weighted_Eigen> matrix_d(
        std::string const &spec_name,
        std::array<SparseSetT *,2> dims,
        RegridParams const &params) const;    // Ignores this->params()

    // ----------- Implements RegridMatrices
    /** Produces its own dims, rather than re-using ones supplied by the user */
    std::unique_ptr<ibmisc::linear::Weighted> matrix(
        std::string const &spec_name) const;

};


}    // namespace
#endif    // guard
