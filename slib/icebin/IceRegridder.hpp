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

#include <unordered_set>
#include <ibmisc/netcdf.hpp>

#include <icebin/AbbrGrid.hpp>
#include <icebin/eigen_types.hpp>
#include <icebin/RegridMatrices.hpp>

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


// ------------------------------------------------------------
class IceRegridder;



/** Represents a single ice sheet.
Produces low-level unscaled matrices *for that single ice sheet*. */
class IceRegridder {
    friend class IceCoupler;
    friend class GCMRegridder_Standard;
    friend class IceWriter;
public:
    typedef GridParameterization Type;


    /** Parent pointer; holds the IceRegridder for ALL ice sheets */
    GCMRegridder const *gcm;
    blitz::Array<double,1> gridA_proj_area;    // Area of GCM's grid cells projected using gridI->sproj

protected:

    Type type;          /// GridParameterization
    std::string _name;  /// "greenland", "antarctica", etc.

public:
    AbbrGrid agridI;            /// Ice grid outlines
    AbbrGrid aexgrid;       /// Exchange grid outlines (between GCM and Ice)
    InterpStyle interp_style;   /// How we interpolate I<-E.  Determines basis functions in E

    // ---------------------------------

    // MatrixFunctions used by corresponding functions in GCMRegridder
    /** Remove unnecessary GCM grid cells. */
    void filter_cellsA(std::function<bool(long)> const &keepA);

public:
    std::string const &name() const { return _name; }

    IceRegridder();
    virtual ~IceRegridder();

//    void clear();

    /**
    @param elevI The elevation of each (unmasked) ice grid cell.
    Indicate a masked-out grid cell with elevI[i] = NaN
    */
    void init(
        std::string const &_name,
        AbbrGrid const &agridA,
        Grid const &_fgridA,
        Grid const &_fgridI,
        Grid const &_fexgrid,
        InterpStyle _interp_style);

    /** @param elevI Elevation at points of ice sheet and land.
    @param maskI Land surface type of each grid cell (ocean, land, ice) */
    void set_elevI(
        DenseArrayT<1> const &_elevI,
        blitz::Array<char,1> const &maskI);

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
    void sApvA(MakeDenseEigenT::AccumT &w) const;

    /** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
    NOTE: wAvAp == sApvA */
    void sEpvE(MakeDenseEigenT::AccumT &w) const;

    /** Produces the unscaled matrix [Interpolation or Ice] <-- [Projected Elevation] */
    virtual void GvEp(MakeDenseEigenT::AccumT &ret,
        blitz::Array<double,1> const *elevI) const = 0;

    /** Produces the unscaled matrix [Interpolation or Ice] <-- [Ice] */
    virtual void GvI(MakeDenseEigenT::AccumT &ret,
        blitz::Array<double,1> const *elevI) const = 0;

    /** Produces the unscaled matrix [Interpolation or Ice] <-- [Projected Atmosphere] */
    virtual void GvAp(MakeDenseEigenT::AccumT &ret,
        blitz::Array<double,1> const *elevI) const = 0;

    /** Define, read or write this data structure inside a NetCDF file.
    @param vname: Variable name (or prefix) to define/read/write it under. */
    virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname);

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
