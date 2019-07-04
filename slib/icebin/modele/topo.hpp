#ifndef ICEBIN_MODELE_TOPO_HPP
#define ICEBIN_MODELE_TOPO_HPP

#include <ibmisc/indexing.hpp>
#include <ibmisc/linear/eigen.hpp>
#include <icebin/modele/hntr.hpp>

/** Utilities involved in computing TOPO files:

   * TOPOO files based on base TOPOO, plus merged EOvAO matrix
   * TOPOA files, based on merged TOPOO file.

NOTES:
 1. GCMRegridder data type is NOT required.

 2. No EXPLICIT merging is done here; the assumption is that
    pre-merged EOvAO matrices alre already provided.

 3. These utilities are used by command-line programs, as well as
    GCMRegridder_ModeleE
*/

namespace icebin {
namespace modele {

/** Encodes relationship between number of ECs in ice model, and
number of ECs in GCM.  GCM has an extra for the "sealand" EC. */
inline int get_nhc_gcm(int nhc_ice)
    { return nhc_ice + 1; }

extern EigenColVectorT compute_wAOm(
blitz::Array<double,1> const &foceanAOp,    // gcmA->foceanAOp
blitz::Array<double,1> const &foceanAOm,    // gcmA->foceanAOm
blitz::Array<double,1> const &wAOp,
SparseSetT &dimAOp,   // Constant
SparseSetT &dimAOm);   // write here


/** Computes EOvEA (raw, unscaled).

Strategy: Group and sum EO elevation classes to produce EA elevation
   classes.  This is done by calling Hntr to discover which grid cells
   in AO correspond to each grid cell in AA.  Weights provided by Hntr
   are ignored, summing up wEO instead.

NOTES:

   1. This strategy assumes that every cell in EO is fully contained
      in a cell of AO.

   2. We should have:
               sum(EOvEA, 1) == wEA

      This can and should be checked; and if it does not work out, we
      should correct by multiplying by:
          EOvEA_corrected = EOvEA_raw * diag(sum(EOvEA,1,'-') .* wEA)
*/
extern void raw_EOvEA(
MakeDenseEigenT::AccumT &&ret,        // {dimEA, dimEO}; dimEO should not change here.
HntrSpec const &hntrO,
HntrSpec const &hntrA,
double const eq_rad,
SparseSetT const *dimAO,            // Used to clip in Hntr::matrix()
blitz::Array<double,1> &wEO_d,            // == EOvI.wM.  Dense indexing.
// Things obtained from gcmA
unsigned int const nhc,    // gcmA->nhc()
ibmisc::Indexing const indexingHCO,    // gcmA->gcmO->indexingHC
ibmisc::Indexing const indexingHCA);    // gcmA->indexingHC



/** Computes EOmvAOm, based on EOpvAOp.  EOpvAOp is converted to EOmvAOm by
removing elements that refer to non-existent grid cells in AOm.
This ultimately provides us with dimEOm as well. */
extern EigenSparseMatrixT compute_EOmvAOm_unscaled(
SparseSetT &dimEOm,        // NOT const
SparseSetT &dimAOm,        // const; pre-computed, should not change
EigenSparseMatrixT const &EOpvAOp,    // Unscaled
SparseSetT const &dimEOp,
SparseSetT const &dimAOp);


extern std::unique_ptr<ibmisc::linear::Weighted_Eigen> _compute_AAmvEAm(
std::array<SparseSetT *,2> dims,
bool scale,        // paramsA.scale
double const eq_rad,    // Radius of the earth

// Things obtained from gcmA
HntrSpec const &hntrO,        // cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr
HntrSpec const &hntrA,        // cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr
ibmisc::Indexing const indexingHCO,    // gcmA->gcmO->indexingHC
ibmisc::Indexing const indexingHCA,    // gcmA->indexingHC
blitz::Array<double,1> const &foceanAOp,    // gcmA->foceanAOp
blitz::Array<double,1> const &foceanAOm,    // gcmA->foceanAOm

// Sub-parts of the computation, pre-computed
EigenSparseMatrixT const &EOpvAOp,
SparseSetT &dimEOp,
SparseSetT &dimAOp,
blitz::Array<double,1> const &wAOp);

/**
@return A textual list of errors found in the sanity check.
*/
std::vector<std::string> make_topoA(
// AAmvEAM is either read from output of global_ec (for just global ice);
// or it's the output of compute_AAmvEAm_merged (for merged global+local ice)
blitz::Array<double,2> const &foceanOp2,
blitz::Array<double,2> const &foceanOm2,     // Rounded FOCEAN
blitz::Array<double,2> const &flakeOm2,
blitz::Array<double,2> const &fgrndOm2,
blitz::Array<double,2> const &fgiceOm2,
blitz::Array<double,2> const &zatmoOm2,
blitz::Array<double,2> const &zlakeOm2,
blitz::Array<double,2> const &zicetopOm2,
//
// Things obtained from gcmA
HntrSpec const &hspecO,        // cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr
HntrSpec const &hspecA,        // cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr
ibmisc::Indexing const indexingHCO,    // gcmA->gcmO->indexingHC   (must reflect local + global ECs)
ibmisc::Indexing const indexingHCA,    // gcmA->indexingHC
std::vector<double> const &hcdefs,        // gcmA->hcdefs()
std::vector<uint16_t> const &underice,    // gcmA->underice
//
double const eq_rad,
EigenSparseMatrixT const &EOpvAOp,        // UNSCALED
SparseSetT &dimEOp,    // const
SparseSetT &dimAOp,    // const
//
//SparseSetT const &dimO,    // Tells us which grid cells in O were changed.
blitz::Array<double,2> &foceanA2,    // Rounded FOCEAN
blitz::Array<double,2> &flakeA2,
blitz::Array<double,2> &fgrndA2,
blitz::Array<double,2> &fgiceA2,
blitz::Array<double,2> &zatmoA2,
blitz::Array<double,2> &zlakeA2,
blitz::Array<double,2> &zicetopA2,
//
blitz::Array<double,3> &fhc3,
blitz::Array<double,3> &elevE3,
blitz::Array<int16_t,3> &underice3);

/** Check that FHC sums to 1 in every gridcell.
@param errors Return any errors found by adding to this vector. */
void sanity_check_fhc(
blitz::Array<double,3> const &fhc3,
std::vector<std::string> &errors);


/** Check that land fractions sum to 1 in every gridcell.
@param errors Return any errors found by adding to this vector. */
void sanity_check_land_fractions(
blitz::Array<double,2> const &foceanA2,    // Rounded FOCEAN
blitz::Array<double,2> const &flakeA2,
blitz::Array<double,2> const &fgrndA2,
blitz::Array<double,2> const &fgiceA2,
std::vector<std::string> &errors);


}}    // namespace
#endif    // guard
