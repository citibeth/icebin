#ifndef ICEBIN_MODELE_MERGE_TOPO_HPP
#define ICEBIN_MODELE_MERGE_TOPO_HPP

#include <ibmisc/blitz.hpp>
#include <ibmisc/zarray.hpp>
#include <ibmisc/indexing.hpp>
#include <icebin/RegridMatrices.hpp>
#include <icebin/GCMRegridder.hpp>

/** Subroutines for merging local ice sheets into base set up based on
global ice (eg from ETOPO1). */

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
@return Any errors encountered on the sanity check
*/
extern void merge_topoO(
// ------ TOPOO arrays, originally with just global ice
// Ice sheets will be merged into them.
//
// Ice model viewpoint (fractional ocean cells)
blitz::Array<double,2> &foceanOp2,    // Fractional FOCEAN
blitz::Array<double,2> &fgiceOp2,
blitz::Array<double,2> &zatmoOp2,
// ModelE viewpoint (rounded ocean cells)
blitz::Array<double,2> &foceanOm2,     // Rounded FOCEAN
blitz::Array<double,2> &flakeOm2,
blitz::Array<double,2> &fgrndOm2,
blitz::Array<double,2> &fgiceOm2,
blitz::Array<double,2> &zatmoOm2,
// Not affected by Om; as long as top of ice is maintained even for ocean-rounded cells.
blitz::Array<double,2> &zicetopO2,
// ------ Local ice to merge in...
GCMRegridder *gcmO,
RegridParams const &paramsA,
std::vector<blitz::Array<double,1>> const &emI_lands,
std::vector<blitz::Array<double,1>> const &emI_ices,
double const eq_rad,    // Radius of the earth
std::vector<std::string> &errors);

/** Return type for compute_EOpvAOp_merged(), ec */
struct EOpvAOpResult {
    SparseSetT dimEOp;    // dimEOp is set and returned; dimAOp is appended
    std::unique_ptr<EigenSparseMatrixT> EOpvAOp;
    std::vector<double> hcdefs; // OUT:  Elev class definitions for merged ice
    ibmisc::Indexing indexingHC;
    std::vector<uint16_t> underice_hc;
};


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
@param hcdefs_base Levels of elevation classes [m] in base EOpvAOp
@param hcdefs OUT: ELevation levels for merged EOpvAOp
@param underice_hc OUT: Expected model underneath the ice of each elevation class
    UI_NOTHING for base ice
    UI_ICEBIN for local ice
@param paramsA Regrid parmaeters that are used by command line programs (known to work)
@return Merged EOpvAOp matrix.  Dense indexing, using dimension maps from dims.
*/
EOpvAOpResult compute_EOpvAOp_merged(  // (generates in dense indexing)
SparseSetT &dimAOp,    // dimAOp is appended
ibmisc::ZArray<int,double,2> const &EOpvAOp_base,    // from linear::Weighted_Compressed
RegridParams paramsO,
GCMRegridder const *gcmO,     // A bunch of local ice sheets
double const eq_rad,    // Radius of the earth
std::vector<blitz::Array<double,1>> const &emI_ices,
bool use_global_ice,
bool use_local_ice,
std::vector<double> const &hcdefs_base, // [nhc]  Elev class definitions for base ice
ibmisc::Indexing const &indexingHC_base,
bool squash_ecs,    // Should ECs be merged if they are the same elevation?
std::vector<std::string> &errors);


/** Merges repeated ECs */
EOpvAOpResult squash_ECs(
std::array<SparseSetT *,2> dims0,
EigenSparseMatrixT const &EOpvAOp0,
std::vector<double> const &hcdefs0, // Elevation of each EC in EOpvAOp0
ibmisc::Indexing const &indexingHC0);



}} // namespace
#endif // guard
