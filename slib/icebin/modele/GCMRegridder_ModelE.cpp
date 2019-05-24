#include <spsparse/eigen.hpp>
#include <ibmisc/stdio.hpp>
#include <ibmisc/linear/compressed.hpp>
#include <icebin/modele/GCMRegridder_ModelE.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/gridgen/GridGen_LonLat.hpp>

using namespace icebin;
using namespace ibmisc;
using namespace spsparse;
using namespace std::placeholders;

namespace icebin {
namespace modele {

static double const NaN = std::numeric_limits<double>::quiet_NaN();

// ----------------------------------------------------------------
// ----------------------------------------------------------------
/** Diagonal matrix converts wOm <- wOp, weight of the elevation
classes in the two different systems (ModelE (m) vs. IceBin (p)).
produces a SCALED matrix.

NOTES:
 1. Meant to be used with class MakeDenseEigenT.

 2. focean_m is fixed over the course of a run, since ModelE can not
    change its ocean mask.

 3. Matrix will be diagonal in sparse indexing.

 4. This is a SCALED matrix, not unscaled; it gives no information on
    the size ("weight") of the grid cells.  It is intended to be used
    only to convert weight vectors between what ModelE vs. the ice
    model sees.

@param foceanAOm
   Ocean surface fraction array (FOCEAN), as seen by
   ModelE.  On Ocean grid, spase indexing.

@param foceanAOp
   Ocean surface fraction array (FOCEAN), as seen by
   ice model.  On Ocean grid, spase indexing.
@param invert
   By default, this produces the matrix AOmvAOp, which converts
   quantities [X m-2] in AOp to quantities [X m-2] in AOm.  In ordert
   to convert cell norms (weights) from AOm<-AOp, one should use
   invert='-', which produces the inverse matrix of invet='+'.
*/
static void scaled_AOmvAOp(
    MakeDenseEigenT::AccumT &&ret,        // {dimAOm, dimAOp}
    blitz::Array<double,1> const &foceanAOp,    // sparse indexing, 0-based
    blitz::Array<double,1> const &foceanAOm,    // sparse indexing, 0-based
    char invert = '+')
{
//    SparseSetT &dimAOp(*ret.dim(1).sparse_set);
    SparseSetT &dimAOp(*ret.sparse_sets[1]);

    // By looping over only things in IceBin ice sheet, we implicitly only
    // look at ice on ocean grid cells where both ModelE and IceBin have something.

    // Look only at ocean grid cells where BOTH ModelE and IcBin have ice
    for (int iAOp_d=0; iAOp_d < dimAOp.dense_extent(); ++iAOp_d) {
        long const iAO_s = dimAOp.to_sparse(iAOp_d);

        double const fcont_p = 1.0 - foceanAOp(iAO_s);
        double const fcont_m = 1.0 - foceanAOm(iAO_s);    // 0 or 1

        if (fcont_m == 0.0) continue;
        if (fcont_m != 1.0) (*icebin_error)(-1,
            "fcont_m[%ld] = %g (%p; %p), must be 0 or 1", iAO_s, fcont_m, foceanAOm.data(), &foceanAOm(iAO_s));

        if (fcont_p == 0.0) continue;    // Can't extend ice that's not there

        // In Om space, cell only accounts for the fraction covered by continent.
        // When divided by wM (scaling), this will have the effect of magnifying
        // cells only partially covered by ocean (in PISM)
        double const val = (invert == '-' ? fcont_p : 1./fcont_p);
        ret.add({iAO_s, iAO_s}, val);
    }
}
// ----------------------------------------------------------------
/** Helper function: clip a cell according to its index */
static bool dim_clip_fn(SparseSetT const *dim, long index)
    { return dim->in_sparse(index); }

static blitz::Array<double,1> dim_clip(SparseSetT const &dim)
{
    // Prepare weight vector to clip
    blitz::Array<double,1> wt(dim.sparse_extent());
    wt = 0.0;
    for (int i_d=0; i_d < dim.dense_extent(); ++i_d) {
        int i_s = dim.to_sparse(i_d);

        wt(i_s) = 1.0;
    }

    return wt;
}

// -------------------------------------------------------------
/** Helper class for raw_EOvEA().
@see raw_EOvEA */
class RawEOvEA {
public:
    MakeDenseEigenT::AccumT ret;    // Final place for EOvEA
    GCMRegridder_ModelE const *gcmA;
    blitz::Array<double,1> const &wEO_d;
    // Things obtained from gcmA
    unsigned int const nhc;    // gcmA->nhc()
    IndexSet const indexingHCO;    // gcmA->gcmO->indexingHC
    IndexSet const indexingHCA;    // gcmA->indexingHC

int n=0;

    RawEOvEA(
        MakeDenseEigenT::AccumT &&_ret,
        GCMRegridder_ModelE const *_gcmA,
        blitz::Array<double,1> const &_wEO_d,
        unsigned int const _nhc,    // gcmA->nhc()
        IndexSet const _indexingHCO,    // gcmA->gcmO->indexingHC
        IndexSet const _indexingHCA)    // gcmA->indexingHC
    : ret(std::move(_ret)), wEO_d(_wEO_d), nhc(_nhc),
        indexingHCO(_indexingHCO), indexingHCA(_indexingHCA) {}

    /** Called by Hntr::matrix() (AvO) */
    void add(std::array<int,2> index, double value)
    {
        // Ignore stray overlaps
        if (std::abs(value) < 1e-8) (*icebin_error)(-1,
            "Found a stray overlap; what should we do about it?");

        SparseSetT &dimEO(*ret.sparse_sets[0]);
        long lAO_s = index[0];
        long lAA_s = index[1];

//if (n < 15) printf("    RawEOvEA: (%d %d) = %g   nc=%d\n", lAO_s, lAA_s, value, gcmA->nhc());
        // Iterate through all possible elevation classes for this gridcell pair
        for (int ihc=0; ihc<nhc; ++ihc) {
            long const lEO_s = indexingHCO.tuple_to_index(
                std::array<long,2>{lAO_s,ihc});

            // Obtain wEO, size of the elevation grid cell
            if (!dimEO.in_sparse(lEO_s)) continue;    // wEO==0 here
            int const lEO_d = dimEO.to_dense(lEO_s);
            double const weightEO = wEO_d(lEO_d);
//if (n < 15) printf("      ihc=%d, lEO_s = %ld   lEO_d=%d  weightEO=%g\n", ihc, lEO_s, lEO_d, weightEO);

            if (weightEO != 0) {
                int const lEA_s = indexingHCA.tuple_to_index(
                    std::array<long,2>{lAA_s,ihc});
                ret.add({lEO_s,lEA_s}, weightEO);
            }
        }
++n;
    }
};

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
static void raw_EOvEA(
    MakeDenseEigenT::AccumT &&ret,        // {dimEA, dimEO}; dimEO should not change here.
    HntrSpec const &hntrO,
    HntrSpec const &hntrA,
    double const eq_rad,
    SparseSetT const *dimAO,            // Used to clip in Hntr::matrix()
    blitz::Array<double,1> &wEO_d,            // == EOvI.wM.  Dense indexing.
    // Things obtained from gcmA
    unsigned int const nhc,    // gcmA->nhc()
    IndexSet const indexingHCO,    // gcmA->gcmO->indexingHC
    IndexSet const indexingHCA)    // gcmA->indexingHC
{
    // Call Hntr to generate AOvAA; and use that (above) to produce EOvEA
    Hntr hntr_AOvAA(17.17, hntrO, hntrA, 0);    // dimB=A,  dimA=O

    hntr_AOvAA.overlap<RawEOvEA, DimClip>(
        RawEOvEA(std::move(ret), wEO_d, nhc, indexingHCO, indexingHCA),
        eq_rad, DimClip(dimAO));
}


// ----------------------------------------------------------------
// ========================================================================

GridSpec_LonLat const &cast_GridSpec_LonLat(GridSpec const &_specO)
{
    // -------- Check types on specO
    GridSpec_LonLat const *specO(dynamic_cast<GridSpec_LonLat const *>(&_specO));
    if (!specO) (*icebin_error)(-1,
        "make_gridA() requires type Grid_LonLat");

    if (specO->north_pole != specO->south_pole) (*icebin_error)(-1,
        "north_pole=%d and south_pole=%d must match in specO",
        specO->north_pole, specO->south_pole);

    if (!specO->hntr.is_set()) (*icebin_error)(-1,
        "make_gridA() requires specO have a Hntr source");

    return *specO;
}




/** Creates Atmosphere grid from an existing Hntr-defined Ocean grid.
It does this creating a Hntr spec for gridA, and then running grid-generation
code on it.

@gridO the Hntr-defined Ocean grid.  gridO->hntr must be set; this
    will happen if it was defined using GridGen_Hntr (or loaded from
    a NetCDF file generated that way).
*/
static AbbrGrid make_agridA(
    std::string const &name,
    SparseSet<long,int> const &dimO,
    HntrSpec const &hspecO,
    double eq_rad)
{

    HntrSpec const hspecA(make_hntrA(hspecO));

    // Use hntr to figure out which grid cells should be realized in A, based
    // on realized grid cells in O
    SparseSet<long,int> dimA;
    Hntr hntrOvA(17.17, hspecO, hspecA, 0);
    HntrGrid const &hgridA(hntrOvA.Bgrid);
    AbbrGrid agridA;
    hntrOvA.overlap(
            accum::SparseSetAccum<SparseSetT,double,2>({nullptr, &dimA}),
        1.0, DimClip(&dimO));

    GridSpec_LonLat specA(make_grid_spec(hspecA, false, 1, eq_rad));
    return make_abbr_grid(name, specA, std::move(dimA));
}

inline AbbrGrid make_agridA(AbbrGrid &agridO)
{
    GridSpec_LonLat const &specO(cast_GridSpec_LonLat(*agridO.spec));
    return make_agridA(agridO.name, agridO.dim, specO.hntr, specO.eq_rad);
}

// ========================================================================
// --------------------------------------------------------
/** Top-level regrid generator computes AOmvAAm.
This is not a regridding matrix ModelE needs directly (it is a component of such).
It exists here as a top-level subroutine so we can test it easily from Python.
@param transpose 'T' if this matrix should be transposed, '.' if not. */
static std::unique_ptr<linear::Weighted_Eigen> compute_AOmvAAm(
    std::array<SparseSetT *,2> dims,
    RegridParams const &paramsA,
    GCMRegridder_ModelE const *gcmA,
    char transpose)
{
    auto &dimAOm(*dims[0]);
    auto &dimAAm(*dims[1]);


    // Obtain hntrA and hntrO for process below.
    GridSpec_LonLat const &specO(cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec));
    GridSpec_LonLat const &specA(cast_GridSpec_LonLat(*gcmA->agridA.spec));

    Hntr hntr_XOmvXAm(17.17, specO.hntr, specA.hntr);

    std::unique_ptr<linear::Weighted_Eigen> ret(new linear::Weighted_Eigen(dims, true));    // conservative
    reset_ptr(ret->M, MakeDenseEigenT(
        std::bind(&Hntr::overlap<MakeDenseEigenT::AccumT,DimClip>,
            &hntr_XOmvXAm, _1, specA.eq_rad, DimClip(&dimAOm)),
        {SparsifyTransform::TO_DENSE_IGNORE_MISSING, SparsifyTransform::ADD_DENSE},
        {&dimAOm, &dimAAm}, transpose).to_eigen());

    // Ignore the scale parameter

    ret->wM.reference(sum(*ret->M, 0, '+'));
    ret->Mw.reference(sum(*ret->M, 1, '+'));

    return ret;
}

// ========================================================================
EigenColVectorT compute_wAOm(
    blitz::Array<double,1> const &foceanAOp,    // gcmA->foceanAOp
    blitz::Array<double,1> const &foceanAOm,    // gcmA->foceanAOm
    blitz::Array<double,1> const &wAOp,
    SparseSetT &dimAOp,   // Constant
    SparseSetT &dimAOm)   // write here
{
    ConstUniverse const_dimAOp({"dimAOp"}, {&dimAOp});

    // ----------- Compute dimAOm properly
    // dimAOm will be a subset of dimAOp.  Remove points in dimAOp that are ocean.
    for (int iAOp_d=0; iAOp_d < dimAOp.dense_extent(); ++iAOp_d) {
        auto const iAO_s(dimAOp.to_sparse(iAOp_d));
        if (foceanAOm(iAO_s) == 0) dimAOm.add_dense(iAO_s);
    }

    // OmvOp
    EigenSparseMatrixT sc_AOmvAOp(MakeDenseEigenT(
        std::bind(&scaled_AOmvAOp, _1, foceanAOp, foceanAOm, '+'),
        {SparsifyTransform::TO_DENSE_IGNORE_MISSING},
        {&dimAOm, &dimAOp}, '.').to_eigen());

    return sc_AOmvAOp * map_eigen_colvector(wAOp);
}
// ========================================================================
/** Computes EOmvAOm, based on EOpvAOp.  EOpvAOp is converted to EOmvAOm by
removing elements that refer to non-existent grid cells in AOm.
This ultimately provides us with dimEOm as well. */
EigenSparseMatrixT compute_EOmvAOm_unscaled(
    SparseSetT &dimEOm,        // NOT const
    SparseSetT &dimAOm,        // const; pre-computed, should not change
    EigenSparseMatrixT const &EOpvAOp,    // Unscaled
    SparseSetT const &dimEOp,
    SparseSetT const &dimAOp)
{
    ConstUniverse const_dimAOm({"dimAOm"}, {&dimAOm});

    // Convert to EOmvAOm by removing cells not in dimAOm
    TupleListT<2> EOmvAOm_tl;
    for (auto ii(begin(EOpvAOp)); ii != end(EOpvAOp); ++ii) {
        auto const iAOp_d = ii->index(1);
        auto const iAOp_s = dimAOp.to_sparse(iAOp_d);
        int iAOm_d;
        if (dimAOm.to_dense_ignore_missing(iAOp_s, iAOm_d)) {
            auto const iEOp_d = ii->index(0);
            auto const iEO_s = dimEOp.to_sparse(iEOp_d);
            auto const iEOm_d = dimEOm.add_dense(iEO_s);  // add to dimEOp
            EOmvAOm_tl.add({iEOm_d, iAOm_d}, ii->value()); // add to EOmvAOm
        }
    }

    // Create the matrix from the TupleList
    EOmvAOm_tl.set_shape(std::array<long,2>{});

    EigenSparseMatrixT EOmvAOm(dimEOm.dense_extent(), dimAOm.dense_extent());
    EOmvAOm.setFromTriplets(EOmvAOm_tl.begin(), EOmvAOm_tl.end());
    return EOmvAOm;
}

// ------------------------------------------------------
/** Computes AAmvEAM (with scaling)
@param dims {dimAAm, dimEAm} User-supplied dimension maps, which will be added to.
@param paramsA Regridding parameters.  NOTE: correctA is ignored.
@param gcmA Used to obtain access to foceanAOp, foceanAOm, gridA, gridO
@param eq_rad Radius of the earth.  Must be the same as eq_rad in ModelE.
@param rmO Regrid matrix generator that can regrid between (AOp, EOp, Ip).
       NOTE: rmO knows these grids as (A, E, I).
*/
static std::unique_ptr<linear::Weighted_Eigen> compute_AAmvEAm_rmO(
    std::array<SparseSetT *,2> dims,
    RegridParams const &paramsA,
    GCMRegridder_ModelE const *gcmA,
    double const eq_rad,    // Radius of the earth
    RegridMatrices_Dynamic const *rmO)
{
    // ------------ Params for generating sub-matrices
    RegridParams paramsO(paramsA);
    paramsO.scale = false;
    paramsO.correctA = true;

    SparseSetT dimEOp, dimAOp;
    std::unique_ptr<linear::Weighted_Eigen> EOpvAOp(
        rmO->matrix_d("EvA", {&dimEOp, &dimAOp}, paramsO));

    return _compute_AAmvEAm(dims, paramsA.scale, eq_rad,
        cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr,
        cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr,
        gcmA->gcmO->indexingHC,
        gcmA->indexingHC,
        gcmA->foceanAOp,
        gcmA->foceanAOm,
        *EOpvAOp->M, dimEOp, dimAOp, EOpvAOp->Mw);
}

// -----------------------------------------------------------
/** Computes some intermediate values used by more than one top-level
    regrid generator */
struct ComputeXAmvIp_Helper {

    /** Allocated stuff that must remain allocated for the duration of
    the temporary values.  This should be transferred to a
    RegridMatrices_Dynamic object being returned by the top-level regrid generator. */
    TmpAlloc tmp;

    // Intermediate regridded values
    SparseSetT dimEOm, dimEOp;

    // Intermediate values
    SparseSetT dimAOp;
    SparseSetT dimAOm;

    std::unique_ptr<EigenSparseMatrixT> XAmvXOm;    // Unscaled matrix

    EigenColVectorT *wXOm_e;
    linear::Weighted_Eigen *XOpvIp;    // TODO: Make this unique_ptr and forgoe the tmp.take() below.
    EigenColVectorT *wXAm_e;

    blitz::Array<double,1> XAmvXOms;

    /** Parameters come from top-level regridding subroutine */
    ComputeXAmvIp_Helper(
        std::array<SparseSetT *,2> dims,
        RegridParams const &paramsA,
        GCMRegridder_ModelE const *gcmA,
        char X,    // Either 'A' or 'E', determines the destination matrix.
        double const eq_rad,    // Radius of the earth
        RegridMatrices_Dynamic const *rmO);
};    // class ComputeXAmvIp_Helper

ComputeXAmvIp_Helper::ComputeXAmvIp_Helper(
    std::array<SparseSetT *,2> dims,
    RegridParams const &paramsA,
    GCMRegridder_ModelE const *gcmA,
    char X,    // Either 'A' or 'E', determines the destination matrix.
    double const eq_rad,    // Radius of the earth
    RegridMatrices_Dynamic const *rmO)
{
    SparseSetT &dimXAm(*dims[0]);
    SparseSetT &dimIp(*dims[1]);

    // ------------ Params for generating sub-matrices
    RegridParams paramsO(paramsA);
    paramsO.scale = false;
    paramsO.correctA = true;

    // We need AOpvIp for wAOP; used below.
    std::unique_ptr<linear::Weighted_Eigen> AOpvIp_corrected(rmO->matrix_d(
        "AvI", {&dimAOp, &dimIp}, paramsO));
    blitz::Array<double,1> const &wAOp(AOpvIp_corrected->wM);

    EigenColVectorT wAOm_e(compute_wAOm(gcmA->foceanAOp, gcmA->foceanAOm, wAOp, dimAOp, dimAOm));

    HntrSpec const &hntrA(cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr);
    HntrSpec const &hntrO(cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr);

    dimXAm.set_sparse_extent(X == 'A' ? gcmA->nA() : gcmA->nE());
    dimIp.set_sparse_extent(rmO->ice_regridder->nI());

    // Ensure that operations don't change universes.  Sanity check.
    std::unique_ptr<ConstUniverse> const_universe;

    if (X == 'E') {
        RegridParams paramsO(paramsA);
        paramsO.scale = false;
        paramsO.correctA = false;

        // Get the universe for EOp / EOm
        const_universe.reset(new ConstUniverse({"dimIp"}, {&dimIp}));
        linear::Weighted_Eigen &EOpvIp(tmp.take(rmO->matrix_d("EvI", {&dimEOp, &dimIp}, paramsO)));

        XOpvIp = &EOpvIp;

        const_universe.reset();        // Check that dimIp hasn't changed

        // ------------------------------------
        SparseSetT dimEOp, dimAOp;
        std::unique_ptr<linear::Weighted_Eigen> EOpvAOp(
            rmO->matrix_d("EvA", {&dimEOp, &dimAOp}, paramsO));

        EigenSparseMatrixT EOmvAOm(compute_EOmvAOm_unscaled(dimEOm, dimAOm, *EOpvAOp->M, dimEOp, dimAOp));

        // wEOm_e
        auto EOmvAOms(sum(EOmvAOm, 1, '-'));
        tmp.take<EigenColVectorT>(wXOm_e,
            EigenColVectorT(EOmvAOm * map_eigen_diagonal(EOmvAOms) * wAOm_e));


        // ---------------- Compute the main matrix
        // ---------------- XAmvIp = XAmvXOm * XOpvIp
        // Rename variables
        SparseSetT &dimEAm(dimXAm);

        blitz::Array<double,1> wXOm(to_blitz(*wXOm_e));
        reset_ptr(XAmvXOm, MakeDenseEigenT(    // TODO: Call this XAvXO, since it's the same 'm' or 'p'
            std::bind(&raw_EOvEA, _1,
                hntrO, hntrA,
                eq_rad, &dimAOm, wXOm,
                gcmA->nhc(), gcmA->gcmO->indexingHC, gcmA->indexingHC),
            {SparsifyTransform::TO_DENSE_IGNORE_MISSING, SparsifyTransform::ADD_DENSE},
            {&dimEOm, &dimEAm}, 'T').to_eigen());
    } else {    // X == 'A'
        RegridParams paramsO(paramsA);
        paramsO.scale = false;
        paramsO.correctA = false;

        auto &AOpvIp(tmp.take<std::unique_ptr<linear::Weighted_Eigen>>(
            rmO->matrix_d("AvI", {&dimAOp, &dimIp}, paramsO)));

        XOpvIp = &*AOpvIp;
        wXOm_e = &wAOm_e;

        // ---------------- Compute the main matrix
        // ---------------- XAmvIp = XAmvXOm * XOpvIp

        SparseSetT &dimAAm(dimXAm);

        // Actually AOmvAAm
        Hntr hntr_XOmvXAm(17.17, hntrO, hntrA);
        reset_ptr(XAmvXOm, MakeDenseEigenT(
            std::bind(&Hntr::overlap<MakeDenseEigenT::AccumT,DimClip>,
                &hntr_XOmvXAm, _1, eq_rad, DimClip(&dimAOm)),
            {SparsifyTransform::TO_DENSE_IGNORE_MISSING, SparsifyTransform::ADD_DENSE},
            {&dimAOm, &dimAAm}, 'T').to_eigen());
    }

    // ---------- wEAm = [scale] * XAmvXOm * wEOm
    // NOTE: We are scaling for the source grid, not destination grid.
    // This will produce a matrix where A = O1+O2+O3+O4, which is what
    // we want for summing up weights of ocean grid.
    XAmvXOms.reference(sum(*XAmvXOm, 1, '-'));    // Not unusual way we weight/scale here

    // The matrix *XAmvXOm * map_eigen_diagonal(XAmvXOms) should consist only of 1's and 0's

    tmp.take(wXAm_e,
        EigenColVectorT(*XAmvXOm * map_eigen_diagonal(XAmvXOms) * *wXOm_e));

}
// -----------------------------------------------------------
/** Given a matrix with one dimension expressed in dense Ap (or WLOG
Op) space...  re-does that dimension to dense Am space.  Assumes that
dimAm is a subset of dimAp.  Dimenions not found will be droped. */
static EigenSparseMatrixT crop_mvp(
    SparseSetT &dimAm, SparseSetT &dimAp,    // (const)
    int index_Apdim,
    EigenSparseMatrixT const &Ap)
{
    TupleListT<2> Am_tl;
    for (auto ii(begin(Ap)); ii != end(Ap); ++ii) {
        auto index(ii->index());
        auto const iAp_d = index[index_Apdim];
        auto const iA_s = dimAp.to_sparse(iAp_d);
        int iAm_d;
        if (dimAm.to_dense_ignore_missing(iA_s, iAm_d)) {
            index[index_Apdim] = iAm_d;
            Am_tl.add(index, ii->value());
        }
    }

    std::array<int,2> shape {Ap.rows(), Ap.cols()};
    shape[index_Apdim] = dimAm.dense_extent();
    EigenSparseMatrixT Am(shape[0], shape[1]);
    Am.setFromTriplets(Am_tl.begin(), Am_tl.end());
    return Am;
}
// -----------------------------------------------------------
/** Computes XAmvIp (with scaling)
@param dims {dimXAm, dimIp} User-supplied dimension maps, which will be added to.
@param paramsA Regridding parameters.  NOTE: correctA is ignored.
@param gcmA Used to obtain access to foceanAOp, foceanAOm, gridA, gridO
@param X Either 'A' or 'E', determines the destination matrix.
@param eq_rad Radius of the earth.  Must be the same as eq_rad in ModelE.
@param rmO Regrid matrix generator that can regrid between (AOp, EOp, Ip).
       NOTE: rmO knows these grids as (A, E, I).
*/
static std::unique_ptr<linear::Weighted_Eigen> compute_XAmvIp(
    std::array<SparseSetT *,2> dims,
    RegridParams const &paramsA,
    GCMRegridder_ModelE const *gcmA,
    char X,    // Either 'A' or 'E', determines the destination matrix.
    double const eq_rad,    // Radius of the earth
    RegridMatrices_Dynamic const *rmO)
{
    std::string dimXA(strprintf("dim%cA", X));
    std::unique_ptr<linear::Weighted_Eigen> ret(new linear::Weighted_Eigen(dims, false));    // not conservative

    // Do the legwork, and fish out things from the results that we need
    ComputeXAmvIp_Helper hh(dims, paramsA, gcmA, X, eq_rad, rmO);
    ret->tmp = std::move(hh.tmp);
    auto &XOpvIp(hh.XOpvIp);
    auto &wXAm_e(*hh.wXAm_e);
    auto &XAmvXOm(hh.XAmvXOm);

    auto &dimXOm(X=='E' ? hh.dimEOm : hh.dimAOm);
    auto &dimXOp(X=='E' ? hh.dimEOp : hh.dimAOp);

printf("XOpvIp: %ld\n", (long)(XOpvIp->M->nonZeros()));
printf("XAmvXOm: %ld\n", (long)(XAmvXOm->nonZeros()));

    // ----------- Put it all together (XAmvIp)
    blitz::Array<double,1> sXOpvIp(1. / XOpvIp->wM);
    ret->wM.reference(to_blitz(wXAm_e));    // wAAm_e or wEAm_e
    if (paramsA.scale) {
        auto sXAm(sum(*XAmvXOm, 0, '-'));
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(sXAm) *  *XAmvXOm *
            crop_mvp(dimXOm, dimXOp, 0,
                map_eigen_diagonal(sXOpvIp) * *XOpvIp->M)));
    } else {
        EigenSparseMatrixT M1(map_eigen_diagonal(sXOpvIp) * *XOpvIp->M);
        EigenSparseMatrixT M2(crop_mvp(dimXOm, dimXOp, 0, M1));
        EigenSparseMatrixT M3(*XAmvXOm * M2);
        ret->M.reset(new EigenSparseMatrixT(*XAmvXOm * M2));
    }
    ret->Mw.reference(XOpvIp->Mw);

    return ret;
}
// -----------------------------------------------------------
/** Computes IpvXAm (with scaling)
@param dims {dimIp, dimXAm} User-supplied dimension maps, which will be added to.
@param paramsA Regridding parameters.  NOTE: correctA is ignored.
@param gcmA Used to obtain access to foceanAOp, foceanAOm, gridA, gridO
@param X Either 'A' or 'E', determines the destination matrix.
@param eq_rad Radius of the earth.  Must be the same as eq_rad in ModelE.
@param rmO Regrid matrix generator that can regrid between (AOp, EOp, Ip).
       NOTE: rmO knows these grids as (A, E, I).
*/
static std::unique_ptr<linear::Weighted_Eigen> compute_IpvXAm(
    std::array<SparseSetT *,2> dims,
    RegridParams const &paramsA,
    GCMRegridder_ModelE const *gcmA,
    char X,    // Either 'A' or 'E', determines the destination matrix.
    double const eq_rad,    // Radius of the earth
    RegridMatrices_Dynamic const *rmO)
{
    std::string dimXA(strprintf("dim%cA", X));
    std::unique_ptr<linear::Weighted_Eigen> ret(new linear::Weighted_Eigen(dims, false));    // not conservative

    // Do the legwork, and fish out things from the results that we need
    ComputeXAmvIp_Helper hh({dims[1], dims[0]}, paramsA, gcmA, X, eq_rad, rmO);
    ret->tmp = std::move(hh.tmp);
    auto &XOpvIp(hh.XOpvIp);
    auto &wXAm_e(*hh.wXAm_e);
    auto &XAmvXOm(hh.XAmvXOm);
    auto &dimXOm(X=='E' ? hh.dimEOm : hh.dimAOm);
    auto &dimXOp(X=='E' ? hh.dimEOp : hh.dimAOp);

    std::string const &matname(X=='A' ? "IvA" : "IvE");
    SparseSetT *dimIp(dims[0]);
    RegridParams paramsO(paramsA);
        paramsO.scale = false;
        paramsO.correctA = false;
    auto IpvXOp(
        rmO->matrix_d(matname,
            {dimIp, &dimXOp}, paramsO));

    // ----------- Put it all together (XAmvIp)
    blitz::Array<double,1> sXOpvIp(1. / XOpvIp->wM);

    auto XOmvXAm(XAmvXOm->transpose());
    auto &sXOmvXAm(hh.XAmvXOms);
    blitz::Array<double,1> sIpvXOp(1. / IpvXOp->wM);
    ret->wM.reference(XOpvIp->Mw);
    if (paramsA.scale) {
        ret->M.reset(new EigenSparseMatrixT(
            crop_mvp(dimXOm, dimXOp, 1,
                map_eigen_diagonal(sIpvXOp) * *IpvXOp->M) *
            map_eigen_diagonal(sXOmvXAm) * XOmvXAm
        ));
    } else {
        ret->M.reset(new EigenSparseMatrixT(
            crop_mvp(dimXOm, dimXOp, 1, *IpvXOp->M) *
            map_eigen_diagonal(sXOmvXAm) * XOmvXAm
        ));
    }
    ret->Mw.reference(to_blitz(wXAm_e));    // wAAm_e or wEAm_e

    return ret;
}
// ============================================================
// Methods for GCMRegridder_ModelE

/** @param gcmO A GCMRegridder_Standard that regrids between AOp,EOp,Ip.
        This is typically loaded directly from a NetCDF file. */
GCMRegridder_ModelE::GCMRegridder_ModelE(
    std::string const &_global_ecO,
    std::shared_ptr<icebin::GCMRegridder> const &_gcmO)
    :  global_ecO(_global_ecO), gcmO(_gcmO)
{
    // Initialize superclass member
    agridA = make_agridA(gcmO->agridA);

    _ice_regridders = &gcmO->ice_regridders();

    correctA = gcmO->correctA;    // Not used
    indexingHC = Indexing(
        {"O", "HC"},
        {0L, 0L},
        {agridA.dim.sparse_extent(), gcmO->indexingHC[1].extent},
        {gcmO->indexingHC.indices()[0], gcmO->indexingHC.indices()[1]});
    indexingE = derive_indexingE(agridA.indexing, indexingHC);

    auto const nO = gcmO->nA();
    foceanAOp.reference(blitz::Array<double,1>(nO));
    foceanAOm.reference(blitz::Array<double,1>(nO));

    // Read number of global EC's out of global EC matrix file.
    if (global_ecO == "") {
        global_hcdefs.clear();
    } else {
        NcIO ncio(global_ecO, 'r');
        ncio_vector(ncio, global_hcdefs, true, "hcdefs", "double",
            get_or_add_dims(ncio, {"nhc"}, {_hcdefs.size()} ));
    }

    // Combine local and global hcdefs
    _hcdefs.clear();
    for (double elev : gcmO->hcdefs()) _hcdefs.push_back(elev);
    for (double elev : global_hcdefs) _hcdefs.push_back(elev);
}




std::unique_ptr<RegridMatrices_Dynamic> GCMRegridder_ModelE::regrid_matrices(
    int sheet_index,
    blitz::Array<double,1> const &elevmaskI,
    RegridParams const &params) const
{
    IceRegridder const *regridder = &*ice_regridders()[sheet_index];
    std::unique_ptr<RegridMatrices_Dynamic> rm(new RegridMatrices_Dynamic(regridder, params));
    GridSpec_LonLat const &specO(cast_GridSpec_LonLat(*gcmO->agridA.spec));

    // Per-ice sheet matrix, makes sense either for:
    //    (a) A matrix involving ice sheet (eg AvI,IvA,EvI,IvE)
    //    (b) Generating global EC's for non-twoway-coupled case, where
    //        the IceRegridder IS the global ice cover
    //        (or at least a shard of it).


    // Delicate construction of rmO
    std::unique_ptr<RegridMatrices_Dynamic> _rmO(
        gcmO->regrid_matrices(sheet_index, elevmaskI, params));
    auto &rmO(rm->tmp.take(std::move(*_rmO)));

    rm->add_regrid("EAmvIp", std::bind(&compute_XAmvIp, _1, _2,
        this, 'E', specO.eq_rad, &rmO));
    rm->add_regrid("AAmvIp", std::bind(&compute_XAmvIp, _1, _2,
        this, 'A', specO.eq_rad, &rmO));


    rm->add_regrid("AAmvEAm", std::bind(&compute_AAmvEAm_rmO, _1, _2,
        this, specO.eq_rad, &rmO));


    rm->add_regrid("IpvEAm", std::bind(&compute_IpvXAm, _1, _2,
        this, 'E', specO.eq_rad, &rmO));
    rm->add_regrid("IpvAAm", std::bind(&compute_IpvXAm, _1, _2,
        this, 'A', specO.eq_rad, &rmO));


    // --------------- Simply-named aliases
    rm->add_regrid("EvI", std::bind(&compute_XAmvIp, _1, _2,
        this, 'E', specO.eq_rad, &rmO));
    rm->add_regrid("AvI", std::bind(&compute_XAmvIp, _1, _2,
        this, 'A', specO.eq_rad, &rmO));


    rm->add_regrid("AvE", std::bind(&compute_AAmvEAm_rmO, _1, _2,
        this, specO.eq_rad, &rmO));


    rm->add_regrid("IvE", std::bind(&compute_IpvXAm, _1, _2,
        this, 'E', specO.eq_rad, &rmO));
    rm->add_regrid("IvA", std::bind(&compute_IpvXAm, _1, _2,
        this, 'A', specO.eq_rad, &rmO));

    // For testing
    rm->add_regrid("AOmvAAm", std::bind(&compute_AOmvAAm, _1, _2,
        this, '.'));
    rm->add_regrid("AAmvAOm", std::bind(&compute_AOmvAAm, _1, _2,
        this, 'T'));

    return rm;
}

}}    // namespace
