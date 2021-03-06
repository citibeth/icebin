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

/** Number of grid cells in the Ice (I) or Exchange (X) grids.
@param gridG Either 'I' or 'X' */
size_t IceRegridder_L0::nG(char gridG) const
{
    switch(gridG) {
        case 'I' : return nI();
        case 'X' : return nX();
        default : (*icebin_error)(-1, "Illegal gridG='%c'", gridG);
     }
}
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
    if (i1 >= n) (*icebin_error)(-1,
        "Elevation %g out of bounds (%g, %g)", xx, xpoints[0], xpoints[xpoints.size()-1]);

    int i0 = i1-1;
    indices[0] = i0;
    indices[1] = i1;
    double ratio = (xx - xpoints[i0]) / (xpoints[i1] - xpoints[i0]);
    weights[0] = (1.0 - ratio);
    weights[1] = ratio;
}




/** Builds an interpolation matrix to go from height points to ice/exchange grid.
@param ret Put the regrid matrix here. */
void IceRegridder_L0::GvEp(
    MakeDenseEigenT::AccumT &&ret,
    char gridG,    // Interpolation grid to use for G: 'I' (ice) or 'G' (exchange)
    blitz::Array<double,1> const *_elevmaskI) const
{
printf("BEGIN IceRegridder_L0::GvEp()\n");
    blitz::Array<double,1> const &elevmaskI(*_elevmaskI);

    if (gcm->hcdefs().size() == 0) (*icebin_error)(-1,
        "IceRegridder_L0::GvEp(): hcdefs is zero-length!");

    // ---------------------------------------
    // Handle Z_INTERP or ELEV_CLASS_INTERP

    // Interpolate in the vertical
    for (int id=0; id<aexgrid.dense_extent(); ++id) {
        long const iA = aexgrid.ijk(id,0);        // GCM Atmosphere grid
        long const iI = aexgrid.ijk(id,1);        // Ice Grid
        long const iX = aexgrid.to_sparse(id);    // X=Exchange Grid
        long const iG = (gridG == 'I' ? iI : iX);   // G=Interpolation Grid

        if (!std::isnan(elevmaskI(iI))) {
            // This cell not masked: look up elevation point as usual
            double const elevation = std::max(elevmaskI(iI), 0.0);

            // Interpolate in height points
            switch(interp_style.index()) {
                case InterpStyle::Z_INTERP :
                {
                    long ihps[2];
                    double whps[2];
                    linterp_1d_b(gcm->hcdefs(), elevation, ihps, whps);

                    if (whps[0] != 0) {
                        auto const iE0 = gcm->indexingHC.tuple_to_index<long,2>({iA, ihps[0]});
                        ret.add({iG, iE0},
                            aexgrid.native_area(id) * whps[0]);
                    }

                    if (whps[1] != 0) {
                        auto const iE1 = gcm->indexingHC.tuple_to_index<long,2>({iA, ihps[1]});
                        ret.add({iG, iE1},
                            aexgrid.native_area(id) * whps[1]);
                    }

                } break;
                case InterpStyle::ELEV_CLASS_INTERP : {
                    int ihps0 = nearest_1d(gcm->hcdefs(), elevation);
                    ret.add({iG, gcm->indexingHC.tuple_to_index<long,2>({iA, ihps0})},
                        aexgrid.native_area(id));
                } break;
            }
        }
    }
printf("END IceRegridder_L0::GvEp()\n");
}
// --------------------------------------------------------
void IceRegridder_L0::GvI(
    MakeDenseEigenT::AccumT &&ret,
    char gridG,    // Interpolation grid to use for G: 'I' (ice) or 'X' (exchange)
    blitz::Array<double,1> const *_elevmaskI) const
{
    blitz::Array<double,1> const &elevmaskI(*_elevmaskI);
    if (gridG == 'I') {
        // Ice <- Ice = Indentity Matrix (scaled)
        // But we need this unscaled... so we use the weight of
        // each grid cell.
        for (int iId=0; iId<agridI.dim.dense_extent(); ++iId) {
            long iIs = agridI.dim.to_sparse(iId);

            // Only include I cells that are NOT masked out
            // Depending on elevmaskI, this could be either just ice
            // or ice and dry land
             if (!std::isnan(elevmaskI(iIs)))
                ret.add({iIs, iIs}, agridI.native_area(iId));
        }
    } else {
        // Exchange <- Ice
        for (int id=0; id<aexgrid.dense_extent(); ++id) {
            // cell->i = index in atmosphere grid
            long const iI = aexgrid.ijk(id,1);        // index in ice grid
            long const iX = aexgrid.to_sparse(id);    // index in exchange grid

            // Only include I cells that are NOT masked out
            // Depending on elevmaskI, this could be either just ice
            // or ice and dry land
            if (!std::isnan(elevmaskI(iI)))
                ret.add({iX,iI}, aexgrid.native_area(id));
        }
    }
}
// --------------------------------------------------------
void IceRegridder_L0::GvAp(
    MakeDenseEigenT::AccumT &&ret,
    char gridG,    // Interpolation grid to use for G: 'I' (ice) or 'X' (exchange)
    blitz::Array<double,1> const *_elevmaskI) const
{
printf("BEGIN IceRegridder_L0::GvAp()\n");
    blitz::Array<double,1> const &elevmaskI(*_elevmaskI);
    for (int id=0; id<aexgrid.dense_extent(); ++id) {
        long const iG = (gridG == 'I' ?
            aexgrid.ijk(id,1) : aexgrid.to_sparse(id));
        long const iA = aexgrid.ijk(id,0);
        long const iI = aexgrid.ijk(id,1);

        // Only include I cells that are NOT masked out
        // Depending on elevmaskI, this could be either just ice
        // or ice and dry land
        if (!std::isnan(elevmaskI(iI)))
            if (aexgrid.native_area(id) > 0) {
                ret.add({iG, iA}, aexgrid.native_area(id));
            }
    }
printf("END IceRegridder_L0::GvAp()\n");
}
// --------------------------------------------------------
void IceRegridder_L0::ncio(NcIO &ncio, std::string const &vname)
{
    IceRegridder::ncio(ncio, vname);
    auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
}

}   // namespace icebin
