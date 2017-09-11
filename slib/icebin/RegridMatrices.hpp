#ifndef ICEBIN_REGRID_MATRICES_HPP
#define ICEBIN_REGRID_MATRICES_HPP

#include <functional>
#include <unordered_set>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>

#include <icebin/eigen_types.hpp>

namespcae icebin {

/** Holds the set of "Ur" (original) matrices produced by an
    IceRegridder for a SINGLE ice sheet. */
class RegridMatrices {
public:
    ibmisc::TmpAlloc tmp;
    class Params;
    typedef std::function<std::unique_ptr<WeightedSparse>(
        std::array<SparseSetT *,2> dims, Params const &params)> MatrixFunction;

    /** Parameters controlling the generation of regridding matrices */
    struct Params {
        /** Produce a scaled vs. unscaled matrix (Scaled means divide
        by the weights vector).  Used for all matrices. */
        bool const scale;

        /** Correct for changes in area due to projections?  Used for
        all matrices. */
        bool const correctA;

        /** If non-zero, smooth the resulting matrix, with sigma as the
        scale length of the smoothing.  Used for IvA and IvE. */
        std::array<double,3> const sigma;

        /** Tells if these parameters are asking us to smooth */
        bool smooth() const { return sigma[0] != 0; }

        Params(bool _scale, bool _correctA, std::array<double,3> const &_sigma) :
            scale(_scale), correctA(_correctA), sigma(_sigma) {}
    };

    std::map<std::string, MatrixFunction> regrids;

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
    @return The regrid matrix and weights
    */
    std::unique_ptr<WeightedSparse> matrix(
        std::string const &spec_name,
        std::array<SparseSetT *,2> dims,
        Params const &_params) const;
};
// -----------------------------------------------------------------
}    // namespace
#endif    // guard
