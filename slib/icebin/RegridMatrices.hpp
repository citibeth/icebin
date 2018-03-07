#ifndef ICEBIN_REGRID_MATRICES_HPP
#define ICEBIN_REGRID_MATRICES_HPP

#include <unordered_set>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>

#include <icebin/eigen_types.hpp>

namespace icebin {

class IceRegridder;


// -----------------------------------------------------------
/** Holds the set of "Ur" (original) matrices produced by an
    IceRegridder for a SINGLE ice sheet. */
class RegridMatrices {
public:
    /** ice_regridder from which this was constructed */
    IceRegridder const * const ice_regridder;

    ibmisc::TmpAlloc tmp;    // Stores local vars for different types.  TODO: Maybe re-do this as simple classmember variables.  At least, see where it used (by removing it and running the compiler)

    typedef std::function<std::unique_ptr<ibmisc::lintransform::Weighted_Eigen>(
        std::array<SparseSetT *,2> dims, RegridParams const &params)> MatrixFunction;

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
    std::unique_ptr<ibmisc::lintransform::Weighted_Eigen> matrix_d(
        std::string const &spec_name,
        std::array<SparseSetT *,2> dims,
        RegridParams const &_params) const;

    // Virtual function promises a matrix in sparse indexing

};
// -----------------------------------------------------------------
/** Used by Python extension to dynamically get regrid matrices in sparse coordinates.  Wraps a RegridMatrices... */
class RegridSet_RegridMatrices : public RegridSet {

    RegridMatrices rm;

public:
    RegridSet_RegridMatrices(RegridMatrices &&_rm, RegridParams const &params) : rm(std::move(_rm)), _params(params) {}


    // ----------- Implements RegridSet
    std::unique_ptr<ibmisc::lintransform::Weighted> matrix(
        std::string const &spec_name,
        RegridParams const &_params) const
    {
        TmpAlloc tmp;
        auto &dims(tmp.make<std::array<SparseSetT,2>>());

        std::unique_ptr<ibmisc::lintransform::Weighted_Eigen>
            M(rm.matrix_d(spec_name, {&dims[0], &dims[1]}, _params));

        M->tmp.merge(std::move(tmp));
        return new std::unique_ptr<ibmisc::lintransform::Weighted>(M.release());
    }

};


std::unique_ptr<ibmisc::lintransform::Weighted>
RegridSet_RegridMatrices::matrix(
    std::string const &spec_name,
    RegridParams const &_params) const
{

    return std::unique_ptr<ibmisc::lintransform::Weighted>(...);
}



}    // namespace
#endif    // guard
