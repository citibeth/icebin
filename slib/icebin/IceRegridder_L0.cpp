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
/** Builds an interpolation matrix to go from height points to ice/exchange grid.
@param ret Put the regrid matrix here.
@param elevIh Must be the result of this->elevI_hash() */
void IceRegridder_L0::GvEp(spsparse::SparseTriplets<SparseMatrix> &ret) const
{
    IceExch dest = interp_grid;
    std::unordered_map<long,double> const elevIh(elevI_hash());


    if (gcm->hpdefs.size() == 0) (*icebin_error)(-1,
        "IceRegridder_L0::GvEp(): hpdefs is zero-length!");

    // ---------------------------------------
    // Handle Z_INTERP or ELEV_CLASS_INTERP

    // Interpolate in the vertical
    for (auto cell = exgrid->cells.begin(); cell != exgrid->cells.end(); ++cell) {
        long const iA = cell->i;        // GCM Atmosphere grid
        long const iI = cell->j;        // Ice Grid
        long const iX = cell->index;    // X=Exchange Grid
        long const iG = (dest == IceExch::ICE ? iI : iX);   // G=Interpolation Grid
        auto ii = elevIh.find(iI);
        if (ii != elevIh.end()) {
            double const elev = ii->second;

            // This cell not masked: look up elevation point as usual
            double elevation = std::max(elev, 0.0);

            // Interpolate in height points
            switch(interp_style.index()) {
                case InterpStyle::Z_INTERP : {
                    int ihps[2];
                    double whps[2];
                    linterp_1d(gcm->hpdefs, elevation, ihps, whps);
                    ret.add({iG, gcm->indexingHP.tuple_to_index<2>({iA, ihps[0]})},
                        cell->native_area * whps[0]);
                    ret.add({iG, gcm->indexingHP.tuple_to_index<2>({iA, ihps[1]})},
                        cell->native_area * whps[1]);
                } break;
                case InterpStyle::ELEV_CLASS_INTERP : {
                    int ihps0 = nearest_1d(gcm->hpdefs, elevation);
                    ret.add({iG, gcm->indexingHP.tuple_to_index<2>({iA, ihps0})},
                        cell->native_area);
                } break;
            }
        }
    }
}
// --------------------------------------------------------
void IceRegridder_L0::GvI(
    spsparse::SparseTriplets<SparseMatrix> &ret) const
{
    std::unordered_map<long,double> const elevIh(elevI_hash());

    if (interp_grid == IceExch::ICE) {
        // Ice <- Ice = Indentity Matrix
        for (auto cell=gridI->cells.begin(); cell != gridI->cells.end(); ++cell) {
            long iI = cell->index;

            if (elevIh.find(iI) != elevIh.end())
                ret.add({iI, iI}, cell->native_area);
        }
    } else {
        // Exchange <- Ice
        for (auto cell = exgrid->cells.begin(); cell != exgrid->cells.end(); ++cell) {
            // cell->i = index in atmosphere grid
            long const iI = cell->j;        // index in ice grid
            long const iX = cell->index;    // index in exchange grid

            if (elevIh.find(iI) != elevIh.end()) {
                ret.add({iX,iI}, cell->native_area);
            }
        }
    }
}
// --------------------------------------------------------
void IceRegridder_L0::GvAp(spsparse::SparseTriplets<SparseMatrix> &ret) const
{
    for (auto cell = exgrid->cells.begin(); cell != exgrid->cells.end(); ++cell) {
        int iG = (interp_grid == IceExch::ICE ? cell->j : cell->index);
        int iA = cell->i;
        if (cell->native_area > 0) ret.add({iG, iA}, cell->native_area);
    }
}
// --------------------------------------------------------
void IceRegridder_L0::ncio(NcIO &ncio, std::string const &vname)
{
    IceRegridder::ncio(ncio, vname);
    auto info_v = get_or_add_var(ncio, vname + ".info", netCDF::ncInt64, {});
    get_or_put_att_enum(info_v, ncio.rw, "interp_grid", interp_grid);
}

}   // namespace icebin
