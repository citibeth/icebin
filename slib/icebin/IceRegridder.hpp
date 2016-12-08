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
#include <spsparse/VectorCooArray.hpp>
#include <icebin/Grid.hpp>

namespace icebin {

class GCMRegridder;
class IceCoupler;

/** Controls how we interpolate from elevation class space to the ice grid */
BOOST_ENUM_VALUES( InterpStyle, int,
    (Z_INTERP)          (0)
    (ELEV_CLASS_INTERP) (1)
)

// ================================================
/** The sparse vector and matrix data types we'll use in IceBin. */
typedef spsparse::VectorCooArray<long, double, 2> SparseMatrix;
typedef spsparse::VectorCooArray<long, double, 1> SparseVector;

typedef spsparse::SparseSet<long,int> SparseSet;
typedef blitz::Array<double,2> DenseMatrix;
typedef blitz::Array<double,1> DenseVector;

typedef Eigen::SparseMatrix<SparseMatrix::val_type> EigenSparseMatrix;

/** Return value of a sparse matrix */
struct WeightedSparse {
    SparseSet dims[2];    // Dense-to-sparse mapping for the dimensions
    std::unique_ptr<EigenSparseMatrix> M;
    DenseVector weight;
};



// ------------------------------------------------------------
/** Represents a single ice sheet.
Produces low-level unscaled matrices *for that single ice sheet*. */
class IceRegridder {
    friend class IceCoupler;

public:
    typedef Grid::Parameterization Type;

    friend class GCMRegridder;
    /** Parent pointer; holds the IceRegridder for ALL ice sheets */
    GCMRegridder const *gcm;

    /** Elevation of grid cells in ice grid (I).
    This also implies a mask: cells with std::isnan() are masked out. */
    blitz::Array<double,1> elevI;
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
        blitz::Array<double,1> &elevI);

    virtual ~IceRegridder();

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
    void sApvA(spsparse::SparseTriplets<SparseMatrix> &w, std::function<bool(long)> const &filter_fn);

    /** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
    NOTE: wAvAp == sApvA */
    void sEpvE(spsparse::SparseTriplets<SparseMatrix> &w, std::function<bool(long)> const &filter_fn);

    /** Produces the unscaled matrix [Interpolation or Ice] <-- [Projected Elevation] */
    virtual void GvEp(
        spsparse::SparseTriplets<SparseMatrix> &ret) const = 0;

    /** Produces the unscaled matrix [Interpolation or Ice] <-- [Ice] */
    virtual void GvI(
        spsparse::SparseTriplets<SparseMatrix> &ret) const = 0;

    /** Produces the unscaled matrix [Interpolation or Ice] <-- [Projected Atmosphere] */
    virtual void GvAp(
        spsparse::SparseTriplets<SparseMatrix> &ret) const = 0;

    /** Define, read or write this data structure inside a NetCDF file.
    @param vname: Variable name (or prefix) to define/read/write it under. */
    virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname);

};  // class IceRegridder

std::unique_ptr<IceRegridder> new_ice_regridder(IceRegridder::Type type);
std::unique_ptr<IceRegridder> new_ice_regridder(ibmisc::NcIO &ncio, std::string const &vname);
// -----------------------------------------------------------
typedef std::function<std::unique_ptr<WeightedSparse>(bool scale, bool correctA)> RegridFunction;

/** Holds the set of "Ur" (original) matrices produced by an
    IceRegridder for a SINGLE ice sheet. */
class RegridMatrices {
public:
    std::map<std::string, RegridFunction> regrids;

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
        bool scale,
        bool correctA) const
    { return (regrids.at(spec_name))(scale, correctA); }
};


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
