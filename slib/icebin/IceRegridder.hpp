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
typedef Eigen::Matrix<val_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> EigenDenseMatrixT;
typedef Eigen::Matrix<val_type, Eigen::Dynamic, 1, Eigen::ColMajor> EigenColVectorT;
typedef Eigen::Matrix<val_type, 1, Eigen::Dynamic, Eigen::ColMajor> EigenRowVectorT;
// -----------------------------------------

/** Return value of a sparse matrix */
struct WeightedSparse {
    std::array<SparseSetT *,2> dims;            // Dense-to-sparse mapping for the dimensions
    std::unique_ptr<EigenSparseMatrixT> M;    // Dense indexing
    DenseArrayT<1> weight;           // Dense indexing

    WeightedSparse(std::array<SparseSetT *,2> _dims) : dims(_dims) {}
};

// ------------------------------------------------------------
class IceRegridder;

/** Holds the set of "Ur" (original) matrices produced by an
    IceRegridder for a SINGLE ice sheet. */
class RegridMatrices {
public:
    class Params;
    typedef std::function<std::unique_ptr<WeightedSparse>(
        std::array<SparseSetT *,2> dims, Params const &params)> Function;

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
        double const sigma;

        Params(bool _scale, bool _correctA, double _sigma) :
            scale(_scale), correctA(_correctA), sigma(_sigma) {}
    };

    std::map<std::string, Function> regrids;

    /** Use this to construct a RegridMatrices instance:
           GCMRegridder gcm_regridder(...);
           auto rm(RegridMatrices(gcm_regridder.sheet("greenland")));
           // rm.regrid("AvI", scale=true, correctA=true)
           auto AvI(rm.regrid("AvI", true, true));
           AvI.M        // The matrix
           AvI.weight   // The weight vector
    */
    RegridMatrices(IceRegridder *sheet);

    /** Retrieves a final regrid matrix.
    @param spec_name: The matrix to produce.
        Should be "AvI", "IvA", "EvI", "IvE", "AvE" or "EvA".
    @param scale: Produce scaled matrix?
        true  --> [kg m-2]
        false --> [kg]
    @param correctA: Correct for projection error in A or E grids?
    @return The regrid matrix and weights
    */
    std::unique_ptr<WeightedSparse> regrid(
        std::string const &spec_name,
        std::array<SparseSetT *,2> dims,
        Params const &params) const
    {
        auto ret(regrids.at(spec_name)(dims, params));
        return ret;
    }
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

    // Functions used by corresponding functions in GCMRegridder
    /** Remove unnecessary GCM grid cells. */
    void filter_cellsA(std::function<bool(long)> const &keepA);

public:
    std::string const &name() const { return _name; }

    IceRegridder();
    void clear();
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

    /** Produces a smoothing matrix on [Ice] grid
    @param sigma Length scale [m] of smoothing. */
    virtual void smoothI(MakeDenseEigenT::AccumT &ret, double sigma) const
        { (*icebin_error)(-1, "IceRegridder::smoothI() not implemented"); }

    /** Define, read or write this data structure inside a NetCDF file.
    @param vname: Variable name (or prefix) to define/read/write it under. */
    virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname, bool rw_full=true);

protected:
    /** Used by RegridMatrices (top-level API) */
    std::unique_ptr<WeightedSparse> compute_IvAE(
        std::array<SparseSetT *,2> dims,
        RegridMatrices::Params const &params, UrAE const &AE);

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
