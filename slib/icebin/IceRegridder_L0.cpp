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

#include <cstdio>
#include <icebin/GCMRegridder.hpp>
#include <icebin/IceRegridder_L0.hpp>

using namespace ibmisc;

namespace icebin {

static double const nan = std::numeric_limits<double>::quiet_NaN();

// --------------------------------------------------------
/** Does elevation-class-style interpolation on height points.  Assumes
elevation class boundaries midway between height points.
@return Index of point in xpoints[] array that is closes to xx. */
static int nearest_1d(
    std::vector<double> const &xpoints,
    double xx)
{
    int n = xpoints.size();

    // This is the point ABOVE our value.
    // (i0 = i1 - 1, xpoints[i0] < xx <= xpoints[i1])
    // See: http://www.cplusplus.com/reference/algorithm/lower_bound/
    int i1 = lower_bound(xpoints.begin(), xpoints.end(), xx) - xpoints.begin();

    // Convert to point NEAREST ot ours
    if (i1 <= 0) return 0;
    else if (i1 >= n) return n-1;
    else {
        int i0 = i1-1;

        // Distance to i0 vs i1
        double d0 = std::abs(xx - xpoints[i0]);
        double d1 = std::abs(xpoints[i1] - xx);

        if (d0 <= d1) return i0;
        return i1;
    }
}


// --------------------------------------------------------
extern void linterp_1d_b(
    std::vector<double> const &xpoints,
    double xx,
    long *indices, double *weights)  // Size-2 arrays
{
    int n = xpoints.size();

    // This is the point ABOVE our value.
    // (i0 = i1 - 1, xpoints[i0] < xx <= xpoints[i1])
    // See: http://www.cplusplus.com/reference/algorithm/lower_bound/
    int i1 = lower_bound(xpoints.begin(), xpoints.end(), xx) - xpoints.begin();

    if (i1 <= 0) i1 = 1;
    if (i1 >= n) i1 = n-1;

    int i0 = i1-1;
    indices[0] = i0;
    indices[1] = i1;
    double ratio = (xx - xpoints[i0]) / (xpoints[i1] - xpoints[i0]);
    weights[0] = (1.0 - ratio);
    weights[1] = ratio;
}




/** Builds an interpolation matrix to go from height points to ice/exchange grid.
@param ret Put the regrid matrix here.
@param elevIh Must be the result of this->elevI_hash() */
void IceRegridder_L0::GvEp(
    MakeDenseEigenT::AccumT &ret,
    blitz::Array<double,1> const *_elevI) const
{
    blitz::Array<double,1> const &elevI(*_elevI);
    IceExch dest = interp_grid;

    if (gcm->hcdefs.size() == 0) (*icebin_error)(-1,
        "IceRegridder_L0::GvEp(): hcdefs is zero-length!");

    // ---------------------------------------
    // Handle Z_INTERP or ELEV_CLASS_INTERP

    // Interpolate in the vertical
    for (auto cell = exgrid->cells.begin(); cell != exgrid->cells.end(); ++cell) {
        long const iA = cell->i;        // GCM Atmosphere grid
        long const iI = cell->j;        // Ice Grid
        long const iX = cell->index;    // X=Exchange Grid
        long const iG = (dest == IceExch::ICE ? iI : iX);   // G=Interpolation Grid

        if (!std::isnan(elevI(iI))) {
            // This cell not masked: look up elevation point as usual
            double elevation = std::max(elevI(iI), 0.0);

            // Interpolate in height points
            switch(interp_style.index()) {
                case InterpStyle::Z_INTERP :
                {
                    long ihps[2];
                    double whps[2];
                    linterp_1d_b(gcm->hcdefs, elevation, ihps, whps);
                    ret.add({iG, gcm->indexingHC.tuple_to_index<long,2>({iA, ihps[0]})},
                        cell->native_area * whps[0]);
                    ret.add({iG, gcm->indexingHC.tuple_to_index<long,2>({iA, ihps[1]})},
                        cell->native_area * whps[1]);
                } break;
                case InterpStyle::ELEV_CLASS_INTERP : {
                    int ihps0 = nearest_1d(gcm->hcdefs, elevation);
                    ret.add({iG, gcm->indexingHC.tuple_to_index<long,2>({iA, ihps0})},
                        cell->native_area);
                } break;
            }
        }
    }
}
// --------------------------------------------------------
void IceRegridder_L0::GvI(
    MakeDenseEigenT::AccumT &ret,
    blitz::Array<double,1> const *_elevI) const
{
    blitz::Array<double,1> const &elevI(*_elevI);
    if (interp_grid == IceExch::ICE) {
        // Ice <- Ice = Indentity Matrix (scaled)
        // But we need this unscaled... so we use the weight of
        // each grid cell.
        for (auto cell=gridI->cells.begin(); cell != gridI->cells.end(); ++cell) {
            long iI = cell->index;

            // Only include I cells that are NOT masked out
            // Depending on elevI, this could be either just ice
            // or ice and dry land
             if (!std::isnan(elevI(iI)))
                ret.add({iI, iI}, cell->native_area);
        }
    } else {
        // Exchange <- Ice
        for (auto cell = exgrid->cells.begin(); cell != exgrid->cells.end(); ++cell) {
            // cell->i = index in atmosphere grid
            long const iI = cell->j;        // index in ice grid
            long const iX = cell->index;    // index in exchange grid

            // Only include I cells that are NOT masked out
            // Depending on elevI, this could be either just ice
            // or ice and dry land
            if (!std::isnan(elevI(iI)))
                ret.add({iX,iI}, cell->native_area);
        }
    }
}
// --------------------------------------------------------
void IceRegridder_L0::GvAp(
    MakeDenseEigenT::AccumT &ret,
    blitz::Array<double,1> const *_elevI) const
{
    blitz::Array<double,1> const &elevI(*_elevI);
    for (auto cell = exgrid->cells.begin(); cell != exgrid->cells.end(); ++cell) {
        long const iG = (interp_grid == IceExch::ICE ? cell->j : cell->index);
        long const iA = cell->i;
        long const iI = cell->j;

        // Only include I cells that are NOT masked out
        // Depending on elevI, this could be either just ice
        // or ice and dry land
        if (!std::isnan(elevI(iI)))
            if (cell->native_area > 0) {
                ret.add({iG, iA}, cell->native_area);
            }
    }
}
// --------------------------------------------------------
void IceRegridder_L0::ncio(NcIO &ncio, std::string const &vname)
{
    IceRegridder::ncio(ncio, vname);
    auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});
    get_or_put_att_enum(info_v, ncio.rw, "interp_grid", interp_grid);
}

}   // namespace icebin
