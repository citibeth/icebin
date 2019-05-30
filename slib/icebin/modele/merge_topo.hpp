#ifndef ICEBIN_MODELE_MERGE_TOPO_HPP
#define ICEBIN_MODELE_MERGE_TOPO_HPP

#include <ibmisc/blitz.hpp>
#include <icebin/ElevMask.hpp>
#include <icebin/RegridMatrices.hpp>
#include <icebin/GCMRegridder.hpp>

namespace icebin {

class GCMRegridder;

namespace modele {

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
// Not affected by Om; as long as top of ice is maintained even for ocean-rounded cells.
blitz::Array<double,2> &zicetopO2,
// ------ Local ice to merge in...
GCMRegridder *gcmO,
RegridParams const &paramsA,
std::vector<ElevMask<1>> const &elevmasks,    // elevation and cover types for each ice sheet
double const eq_rad);    // Radius of the earth


EigenSparseMatrixT compute_EOpvAOp_merged(  // (generates in dense indexing)
std::array<SparseSetT *,2> dims,
std::string const &global_ecO,    // File written by global_ec
RegridParams const &paramsA,
GCMRegridder const *gcmO,     // A bunch of local ice sheets
double const eq_rad,    // Radius of the earth
std::vector<ElevMask<1>> const elevmasks,    // elevation and cover types for each ice sheet
bool use_global_ice,
bool use_local_ice,
bool include_bedrock);    // true if non-ice covered areas of land should also be included




}} // namespace
#endif // guard
