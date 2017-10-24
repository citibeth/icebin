#ifndef ICEBIN_MODELE_GCMREGRIDDER_MODELE_HPP
#define ICEBIN_MODELE_GCMREGRIDDER_MODELE_HPP

#include <icebin/GCMRegridder.hpp>

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
    blitz::Array<double,1> _foceanAOm, _foceanAOp;    // Hold allocated memory here

public:
    /** A GCMRegridder_Standard that regrids between AOp,EOp,Ip.
    This is typically loaded directly from a NetCDF file. */
    std::shared_ptr<icebin::GCMRegridder> const gcmO;

    /** ModelE ocean cover, on the Ocean grid, as seen by the ice
    model.  Ocean grid cells can contain fractional ocean
    cover.  foceanAOp can change over the course of a ModelE run,
    as the ice model's ice extent changes.
    @see ModelE FOCEAN or FOCEN */
    blitz::Array<double,1> foceanAOp;

    /** ModelE ocean cover, on the Ocean grid, as seen by ModelE.
    Cells are either all ocean (==1.0) or all continent (==0.0).
    foceanAOm does NOT change over the course of a ModelE run,
    because the ModelE ocean is not able to change shape mid-run. */
    blitz::Array<double,1> foceanAOm;

    /** @param _gcmO Underlying regridder for the Ocean Grid Regime. */
    GCMRegridder_ModelE(std::shared_ptr<icebin::GCMRegridder> const &_gcmO);

    /**
    References to these variables must be kept alive externally.
    NOTE: Used only by the Python interface.
    @param _foceanAOm Ocean surface fraction array (FOCEAN), as seen by
       ModelE; sparse indexing.  On Ocean grid.
    @param _foceanAOp Ocean surface fraction array (FOCEAN), as seen by
       ice model; sparse indexing.  On Ocean grid. */
    void set_focean(
        blitz::Array<double,1> &_foceanAOp,
        blitz::Array<double,1> &_foceanAOm);

    // ------------------------------------------------------------
    // Override virtual functions
    IceRegridder *ice_regridder(std::string const &name) const;

    /** Constructs the top-level transformations this regrid generator
    knows how to make.  They are: EAmvIp, AAmvIp, AAmvEAm, IpvEAm, IpvAAm.
    NOTES:
       1. Grid names here are ModelE_specific, not the simple A,E,I generated
          by GCMRegridder_Standard.
       2. It can also generate AOmvAAm and AAmvAOm, for testing purposes. */
    RegridMatrices regrid_matrices(
        std::string const &ice_sheet_name,
        blitz::Array<double,1> const *elevI) const;
};

}}    // namespace
#endif
