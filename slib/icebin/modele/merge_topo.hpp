#ifndef ICEBIN_MODELE_MERGE_TOPO_HPP
#define ICEBIN_MODELE_MERGE_TOPO_HPP

#include <ibmisc/blitz.hpp>
#include <ibmisc/zarray.hpp>
#include <icebin/RegridMatrices.hpp>
#include <icebin/GCMRegridder.hpp>

namespace icebin {

class GCMRegridder;

namespace modele {

/** Merge per-ice sheet data into a global base TOPO files (from which
ice sheet(s) have be removed; output of topo_base.cpp).  Some fields
from input TOPO file need to be in "raw" and "rounded" form.  Rounding
means that ocean grid cells have been rounded to be all-ocean or
no-ocean; and other fields have been adjusted accordingly.  This is a
core requirement of the ModelE ocean.

@param gcmO The GCMRegridder that will provide relevant data for the
    (local) ice sheets to be merged in.
@param paramsA Regrid parameters to use when processing local ice sheets.
@param emI_lands Per-ice sheet elevmaskI arrays (elevation where there
    is ice or land; NaN where there is not.)
@param emI_ices Per-ice sheet elevmaskI arrays (elevation where there
    is ice; NaN where there is not.)
*/
extern void merge_topoO(
// ------ TOPOO arrays, originally with just global ice
// Ice model viewpoint (fractional ocean cells)
blitz::Array<double,2> &foceanOp2,    // Fractional FOCEAN
blitz::Array<double,2> &fgiceOp2,
blitz::Array<double,2> &zatmoOp2,
// ModelE viewpoint (rounded ocean cells)
blitz::Array<double,2> &foceanOm2,     // Rounded FOCEAN
blitz::Array<double,2> &fgrndOm2,
blitz::Array<double,2> &fgiceOm2,
blitz::Array<double,2> &zatmoOm2,
// Not affected by ocean gridcell rounding; as long as top of ice is maintained even for ocean-rounded cells.
blitz::Array<double,2> &zicetopO2,
// ------ Local ice to merge in...
GCMRegridder *gcmO,
RegridParams const &paramsA,
std::vector<blitz::Array<double,1>> const &emI_lands,
std::vector<blitz::Array<double,1>> const &emI_ices,
double const eq_rad);    // Radius of the earth



/** Merge per-ice sheet data into a global base EOpvAOp matrix (from
which ice sheet(s) have be removed; output of global_ec.cpp).  The
matrix is all based on un-rounded ("raw") verions of TOPO fields.

@param dims Dimension maps, to be added to as needed when generating
    merged EOpvAOp matrix.
@param EOpvAOP_base EOpvAOp matrix, based on global ice only.  Sparse
    indexing, so no dimension maps needed.
@param emI_ices Per-ice sheet elevmaskI arrays (elevation where there
    is ice; NaN where there is not.)
@param use_global_ice Merge global ice into output matrix?  Normally
    true.  Set to false in order to generate EOpvAOp matrix from local
    ice only.
@param use_local_ice Merge local ice into output matrix?  Normally true.
@return Merged EOpvAOp matrix.  Dense indexing, using dimension maps from dims.
*/
EigenSparseMatrixT compute_EOpvAOp_merged(  // (generates in dense indexing)
std::array<SparseSetT *,2> dims,
ibmisc::ZArray<int,double,2> const &EOpvAOp_base,    // from linear::Weighted_Compressed
RegridParams const &paramsA,
GCMRegridder const *gcmO,     // A bunch of local ice sheets
double const eq_rad,    // Radius of the earth
std::vector<blitz::Array<double,1>> const &emI_ices,
bool use_global_ice,
bool use_local_ice);

}} // namespace
#endif // guard
