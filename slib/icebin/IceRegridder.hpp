/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ICEBIN_ICEREGRIDDER_H
#define ICEBIN_ICEREGRIDDER_H

#include <functional>
#include <unordered_set>
#include <ibmisc/netcdf.hpp>
#include <spsparse/eigen.hpp>

#include <icebin/Grid.hpp>

namespace icebin {

class GCMRegridder;
class IceCoupler;
class IceWriter;    // Adjoint to IceCoupler

/** Controls how we interpolate from elevation class space to the ice grid */
BOOST_ENUM_VALUES( InterpStyle, int,
    (Z_INTERP)          (0)
    (ELEV_CLASS_INTERP) (1)
)

// ================================================

// Types that will be used throughout as template arguments
typedef long sparse_index_type;
typedef int dense_index_type;
typedef double val_type;

// -----------------------------------------
typedef spsparse::MakeDenseEigen<sparse_index_type, val_type, 0, dense_index_type> MakeDenseEigenT;
template<int RANK>
    using TupleListT = MakeDenseEigenT::TupleListT<RANK>;
template<int RANK>
    using DenseArrayT = blitz::Array<val_type,RANK>;
typedef MakeDenseEigenT::SparseSetT SparseSetT;
typedef MakeDenseEigenT::EigenSparseMatrixT EigenSparseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic> EigenDenseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, 1> EigenColVectorT;
typedef Eigen::Matrix<val_type, 1, Eigen::Dynamic> EigenRowVectorT;
typedef Eigen::DiagonalMatrix<val_type, Eigen::Dynamic> EigenDiagonalMatrixT;
// -----------------------------------------

/** Return value of a sparse matrix */
struct WeightedSparse {
    std::array<SparseSetT *,2> dims;            // Dense-to-sparse mapping for the dimensions


    // If M=BvA, then wM = wBvA = area of B cells
    DenseArrayT<1> wM;           // Dense indexing

    std::unique_ptr<EigenSparseMatrixT> M;    // Dense indexing
    bool smooth;    // Is M smoothed?

    /** Should we enforce conservation (in the face of smoothing)?
    Unsmoothed regrid matrices are already conservative... */
    bool conserve;

    // Area of A cells
    DenseArrayT<1> Mw;

    WeightedSparse(std::array<SparseSetT *,2> _dims) : dims(_dims), smooth(false) {}

    /** Applies a regrid matrix.
    Nominally computes B{in} = smoothB{ii} * BvA{ij} * A{jn}
    (where BvA is this)
    In the face of smoothing, it also compensates for non-conservation in
    smoothB.

        |i| = Size of B vector space
        |j| = Size of A vector space
        |n| = Number of vectors being transformed

    @param A The values to regrid, as a series of Eigen column vectors.
    */
    EigenDenseMatrixT apply_e(
        // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable
        double fill = std::numeric_limits<double>::quiet_NaN()) const;    // Fill value for cells not in BvA matrix


    /** Apply to multiple variables */
    blitz::Array<double,2> apply(
        // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable
        double fill = std::numeric_limits<double>::quiet_NaN(),    // Fill value for cells not in BvA matrix
        TmpAlloc &tmp) const
    {
        return to_blitz(apply_e(A_b, fill), tmp);
    }


    /** Apply to a single variable */
    blitz::Array<double,1> apply(
        // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
        blitz::Array<double,1> const &A_b,       // A_b{j} One variable
        double fill = std::numeric_limits<double>::quiet_NaN(),    // Fill value for cells not in BvA matrix
        TmpAlloc &tmp) const
    {
        auto A_b2(reshape<double,1,2>(A_b, {1, A_b.shape()[0]}));
        auto ret2(to_blitz(apply_e(A_b2, fill), tmp));
        return reshape<double,2,1>(ret2, {ret2.shape()[1]});
    }


    /** Read/write to NetCDF */
    void ncio(ibmisc::NcIO &ncio,
        std::string const &vname,
        std::array<std::string,2> dim_names);

};

// ------------------------------------------------------------
class IceRegridder;

/** Holds the set of "Ur" (original) matrices produced by an
    IceRegridder for a SINGLE ice sheet. */
class RegridMatrices {
public:
    class Params;
    typedef std::function<std::unique_ptr<WeightedSparse>(
        std::array<SparseSetT *,2> dims, Params const &params)> MatrixFunction;

    typedef std::function<void(
        TupleListT<2> &ret,
        SparseSetT const &dimX,
        DenseArrayT<1> const &area_d,    // RM.regrid()->wM
        std::array<double,3> const &sigma
        )> SmoothingFunction;

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

        /** Should we enforce conservation (in the face of smoothing)?
        Unsmoothed regrid matrices are already conservative... */
        bool const conserve;

        bool smooth() const { return sigma[0] != 0; }

        Params(bool _scale, bool _correctA, std::array<double,3> const &_sigma, bool _conserve=false) :
            scale(_scale), correctA(_correctA), sigma(_sigma), conserve(_conserve) {}
    };

protected:
    struct Binding {
        MatrixFunction const matrix;
        SmoothingFunction const *smooth;

        Binding(MatrixFunction const &_matrix, SmoothingFunction const *_smooth) :
            matrix(_matrix), smooth(_smooth) {}
    };

    // Smoothing functions for the different grids
    RegridMatrices::SmoothingFunction smoothI;

    std::map<std::string, Binding> regrids;
    void add_regrid(std::string const &spec,
        SmoothingFunction const *smooth,
        MatrixFunction const &regrid);


public:
    /** Use this to construct a RegridMatrices instance:
           GCMRegridder gcm_regridder(...);
           auto rm(RegridMatrices(gcm_regridder.sheet("greenland")));
           // rm.regrid("AvI", scale=true, correctA=true)
           auto AvI(rm.regrid("AvI", true, true));
           AvI.M        // The matrix
           AvI.wM       // The weight vector
    */
    RegridMatrices(IceRegridder *sheet);

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
class UrAE;

/** Represents a single ice sheet.
Produces low-level unscaled matrices *for that single ice sheet*. */
class IceRegridder {
    friend class IceCoupler;
    friend class RegridMatrices;
public:
    typedef Grid::Parameterization Type;

    friend class GCMRegridder;
    friend class IceWriter;

    /** Parent pointer; holds the IceRegridder for ALL ice sheets */
    GCMRegridder const *gcm;

    /** Elevation of grid cells in ice grid (I).
    This also implies a mask: cells with std::isnan() are masked out. */
    DenseArrayT<1> elevI;
protected:

    Type type;
    std::string _name;  /// "greenland", "antarctica", etc.
    std::unique_ptr<Grid> gridI;            /// Ice grid outlines
    std::unique_ptr<Grid> exgrid;       /// Exchange grid outlines (between GCM and Ice)
    InterpStyle interp_style;   /// How we interpolate I<-E

    // ---------------------------------

    // MatrixFunctions used by corresponding functions in GCMRegridder
    /** Remove unnecessary GCM grid cells. */
    void filter_cellsA(std::function<bool(long)> const &keepA);

public:
    std::string const &name() const { return _name; }

    IceRegridder();
    void clear();

    /**
    @param elevI The elevation of each (unmasked) ice grid cell.
    Indicate a masked-out grid cell with elevI[i] = NaN
    */
    void init(
        std::string const &_name,
        std::unique_ptr<Grid> &&_gridI,
        std::unique_ptr<Grid> &&_exgrid,
        InterpStyle _interp_style,
        DenseArrayT<1> &elevI);

    virtual ~IceRegridder();

    void set_elevI(DenseArrayT<1> const &_elevI)
        { elevI = _elevI; }

    // ------------------------------------------------
    /** Number of dimensions of ice vector space */
    virtual size_t nI() const = 0;

    /** Number of dimensions of interpolation grid vector space. */
    virtual size_t nG() const = 0;

    // ------------------------------------------------
    // Matrix production subroutines here store their output in
    // spsparse accumulators (Anything that accepts a series of
    // (i,j,value) triplets.  The ``SparseTriplets`` class is a
    // special-purpose accumulator that can then produce Eigen
    // matrices as output.


    /** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
    NOTE: wAvAp == sApvA */
    void sApvA(MakeDenseEigenT::AccumT &w);

    /** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
    NOTE: wAvAp == sApvA */
    void sEpvE(MakeDenseEigenT::AccumT &w);

    /** Produces the unscaled matrix [Interpolation or Ice] <-- [Projected Elevation] */
    virtual void GvEp(MakeDenseEigenT::AccumT &ret) const = 0;

    /** Produces the unscaled matrix [Interpolation or Ice] <-- [Ice] */
    virtual void GvI(MakeDenseEigenT::AccumT &ret) const = 0;

    /** Produces the unscaled matrix [Interpolation or Ice] <-- [Projected Atmosphere] */
    virtual void GvAp(MakeDenseEigenT::AccumT &ret) const = 0;

    /** Define, read or write this data structure inside a NetCDF file.
    @param vname: Variable name (or prefix) to define/read/write it under. */
    virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname, bool rw_full=true);
};  // class IceRegridder

std::unique_ptr<IceRegridder> new_ice_regridder(IceRegridder::Type type);
std::unique_ptr<IceRegridder> new_ice_regridder(ibmisc::NcIO &ncio, std::string const &vname);
// -----------------------------------------------------------


/** Gives weights for linear interpolation with a bunch of points.  If
our point is off the end of the range, just continue the slope in
extrapolation.

For example, if:
    xpoints = {3, 5, 6, 8}
    xx = 6.5
Then this subroutine will set:
    indices = {2,3}
    weights = {.75, .25}

@param xpoints The points between which we will interpolate.
    NOTE: this is not blitz::Array<double,1> because Blitz++ does not
          (yet) implement STL-compatible iterators.
@param xx The point for which we want an interpolation formula
@param indices Place to store indices for the interpolation formula.
    Array of exactly two elements.
@param weights Place to store weights for the interpolation formula.
    Array of exactly two elements.
*/
extern void linterp_1d(
    std::vector<double> const &xpoints,
    double xx,
    int *indices, double *weights); // Size-2 arrays

}    // namespace
#endif    // guard
