#ifndef ICEBIN_REGRID_MATRICES_HPP
#define ICEBIN_REGRID_MATRICES_HPP

#include <unordered_set>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>

#include <icebin/eigen_types.hpp>

namespace icebin {

class IceRegridder;


/** Return value of a sparse matrix */
struct WeightedSparse {
    ibmisc::TmpAlloc tmp;            // Stuff that needs same lifetime as WeightedSparse

    std::array<SparseSetT *,2> dims;            // Dense-to-sparse mapping for the dimensions

    // If M=BvA, then wM = wBvA = area of B cells
    DenseArrayT<1> wM;           // Dense indexing

    std::unique_ptr<EigenSparseMatrixT> M;    // Dense indexing

    // Area of A cells
    DenseArrayT<1> Mw;


    /** True if this regridding matrix is conservative.  Matrices could be
    non-conservative, for example, in the face of smoothing on I.  Or when
    regridding between the IceBin and ModelE ice sheets. */
    bool conservative;

    WeightedSparse(std::array<SparseSetT *,2> _dims, bool _conservative) : dims(_dims), conservative(_conservative) {}

    /** Applies a regrid matrix.
    Nominally computes B{in} = BvA{ij} * A{jn}
    (where BvA is this)
    If the matrix ix non-conservative, it also accounts for that.

        |i| = Size of B vector space
        |j| = Size of A vector space
        |n| = Number of vectors being transformed

    @param A The values to regrid, as a series of Eigen column vectors.
    @return Eigen type
    */
    EigenDenseMatrixT apply_e(
        // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable (_b means dense here)
        double fill = std::numeric_limits<double>::quiet_NaN(),    // Fill value for cells not in BvA matrix
        bool force_conservation=true) const;     // Set if you want apply_e() to conserve, even if !M->conservative


    /** Apply to multiple variables
    @return Blitz type */
    blitz::Array<double,2> apply(
        // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable
        double fill,    // Fill value for cells not in BvA matrix
        bool force_conservation,
        ibmisc::TmpAlloc &tmp) const;


    /** Apply to a single variable */
    blitz::Array<double,1> apply(
        // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,1> const &A_b,       // A_b{j} One variable
        double fill,    // Fill value for cells not in BvA matrix
        bool force_conservation,
        ibmisc::TmpAlloc &tmp) const;

    /** Read/write to NetCDF */
    void ncio(ibmisc::NcIO &ncio,
        std::string const &vname,
        std::array<std::string,2> dim_names);

};
// -----------------------------------------------------------
/** Holds the set of "Ur" (original) matrices produced by an
    IceRegridder for a SINGLE ice sheet. */
class RegridMatrices {
public:
    /** ice_regridder from which this was constructed */
    IceRegridder const * const ice_regridder;

    ibmisc::TmpAlloc tmp;    // Stores local vars for different types.  TODO: Maybe re-do this as simple classmember variables.  At least, see where it used (by removing it and running the compiler)

    /** Parameters controlling the generation of regridding matrices */
    struct Params {
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

    typedef std::function<std::unique_ptr<WeightedSparse>(
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
