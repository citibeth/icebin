#ifndef ICEBIN_MODELE_GCMREGRIDDER_MODELE_HPP
#define ICEBIN_MODELE_GCMREGRIDDER_MODELE_HPP

#include <ibmisc/linear/eigen.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/modele/grids.hpp>

namespace icebin {
namespace modele {

/** This class produces regridding matrices between the grids (AAm,
EAm, Ip).

Regridding Matrices
-------------------

These matrices are computed based on another GCMRegridder that
regrids between (AOp, EOp, Ip).  Matrix names are constructed as
follows:

Xp = Grid X, as seen by the dynamic ice model (p = "PISM")
Xm = Grid X, as seen by the GCM (m = "ModelE")

For the Grid named XY:
  * First letter: Grid, as defined by "traditional" GCMRegridder:
      A = Atmosphere Grid
      I = Ice Grid
      E = Elevation Grid (intermediated between A and I)

  * Second letter: Grid Regime; the definition of the "atmosphere"
      grid in a particular GCMRegridder.

      O = Ocean Grid Regime.  The underlying regridder defines
          grids all based on the "atmosphere" grid actually being
          the ModelE ocean grid.  To be useful for the actual
          ModelE, ocean grid cells in matrices from this regime
          must be combined into their corresponding atmosphere
          grid cells (@see raw_AOvAA(), raw_EOvEA()).

      A = Atmosphere Grid Regime.  This GCMRegridder defines grids
          based on the "atmosphere" being the real ModelE
          atmosphere grid.  This is done by using gcmO (The Ocean
          Grid Regime) to produce matrices, and then adjusting
          them for the atmopshere grid.

Putting it all together, we deal with the following grids:

    AAm = ModelE's Atmosphere Grid
    EAm = ModelE's Elevation Grid
    Ip = Dynamic Ice Model Grid
    AOm = ModelE's Ocean Grid
    AOp = Ice Model's Ocean Grid (has different ice extent than in ModelE)
    EOm = Elevation grid based on ModelE's Ocean Grid
    EOp = Elevation grid based on Ice Model's Ocean Grid


Weight Vectors
--------------

A weight vector defines the "weight" (integral) of each basis function
in a vector space.  Weight vectors that go with a regridding matrix
are named according to the matrix itself:

    wBvA = Weight vector for grid B from the matrix BvA
    BvAw = Weight vector for grid A from the matrix BvA

Sometimes, these are abbreviated as wB or wA.

Weight vector variables may sometimes have an _e suffix, if they are
stored in an Eigen-based data structure; un-suffixed weight vectors
are blitz::Array<double,1>.

Indexing
--------

Matrices are computed in a re-numbered subspace of each vectors space
(one that eliminates unused cells).  Index variables into vector
spaces have suffixes that determine whether this index refers to the
original grid ("sparse indexing") or the derived grid ("dense
indexing").  The SparseSetT variables convert between the two kinds of
indexing.

  * Suffix: Sparse vs. Dense indexing

      _s = Sparse indexing
           The natural 1-D indexing for a grid.
           Sparse-indexed vectors may contain many unused elements.

      _d = Dense indexing
           Indices are produced only for grid cells that are used.
           Dense <-> Sparse conversion provided by a SparseSetT.
           Dense-indexed vectors use all their elements.

    Matrices are sparse, whether sparse or dense indexing is used.
    However, multiplying sparse-indexed sparse matrices takes more
    memory than multiplying dense-indexed sparse matrices.

*/
class GCMRegridder_ModelE : public GCMRegridder {
//    blitz::Array<double,1> _foceanAOm, _foceanAOp;    // Hold allocated memory here

public:


    /** Name of file created by global_ec, used for matrices on global
    (non-ice shet) ice with elevation classes.  Generating those
    matrices on-the-fly is not practical because the global ice grid
    is REALLY big.  Instead, we read EvA/AvE out of a pre-computed
    file.  The global_ec file MUST be generated, on the OCEAN grid,
    with the following RegridParams:
        scale = false
        correctA = true

    NOTE: if (global_ecO==""), then we should just use the ice
          sheet-supplied GCMRegridderes (gcmO).  This mode is used by
          global_ec.cpp when computing global ECs for direct use in
          ModelE, without two-way ice sheet coupling.
    */
    std::string const global_ecO;

    std::vector<double> global_hcdefs;

public:
    /** A GCMRegridder_Standard that regrids between AOp,EOp,Ip for ice sheets.
    This is typically loaded directly from a NetCDF file. */
    std::shared_ptr<icebin::GCMRegridder> const gcmO;

    /** Base EOpvAOp matrix, laoded from TOPO_OC file */
    ibmisc::ZArray<int,double,2> EOpvAOp_base,    // from linear::Weighted_Compressed

#if 0
    /** ModelE ocean cover, on the Ocean grid, as seen by the ice
    model (sparse indexing).  Ocean grid cells can contain fractional
    ocean cover.  foceanAOp can change over the course of a ModelE
    run, as the ice model's ice extent changes.
    @see ModelE FOCEAN or FOCEN */
    blitz::Array<double,1> foceanAOp;

    /** ModelE ocean cover, on the Ocean grid, as seen by ModelE
    (sparse indexing).  Cells are either all ocean (==1.0) or all
    continent (==0.0).  foceanAOm does NOT change over the course of a
    ModelE run, because the ModelE ocean is not able to change shape
    mid-run. */
    blitz::Array<double,1> foceanAOm;
#endif

    /** Constructor used in coupler: create the GCMRegridder first,
        then fill in foceanAOp and foceanAOm later.
    @param _gcmO Underlying regridder for the Ocean Grid Regime. */
    GCMRegridder_ModelE(
        std::string const &_global_ecO,
        std::shared_ptr<icebin::GCMRegridder> const &_gcmO);

    HntrSpec const &hspecO()
        { return cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr; }
    HntrSpec const &hspecA()
        { return cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr; }
    GridSpec_LonLat const &specO()
        { return cast_GridSpec_LonLat(*gcmO->agridA.spec); }

TODO: get rid of eq_rad in metaO

    /** Determines whether an elevation class is handled by IceBin or
    ModelE push-down */
    uint16_t underice(int ihc)
        { return (ihc < gcmO->nhc() ? UI_ICEBIN : UI_NOTHING); }

    // ------------------------------------------------------------
    // Override virtual functions
    IceRegridder *ice_regridder(std::string const &name) const;

    /** Constructs the top-level transformations this regrid generator
    knows how to make.  They are: EAmvIp, AAmvIp, AAmvEAm, IpvEAm, IpvAAm.
    NOTES:
       1. Grid names here are ModelE_specific, not the simple A,E,I generated
          by GCMRegridder_Standard.
       2. It can also generate AOmvAAm and AAmvAOm, for testing purposes. */
    std::unique_ptr<RegridMatrices_Dynamic> regrid_matrices(
        int sheet_index,
        blitz::Array<double,1> const &elevmaskI,
        RegridParams const &params = RegridParams()) const;

    linear::Weighted_Tuple global_AvE(
        std::vector<blitz::Array<double,1>> const &emI_lands,
        std::vector<blitz::Array<double,1>> const &emI_ices,
        RegridParams const &params) const;

};

/** Casts to a Grid_Lonlat, which is what we know is used by ModelE */
GridSpec_LonLat const &cast_GridSpec_LonLat(GridSpec const &_specO);


// ====================================================================
// Stuff used by make_topo.cpp





}}    // namespace
#endif
