#include <icebin/modele/GCMRegridder_ModelE.hpp>
#include <icebin/modele/hntr.hpp>

using namespace icebin;
using namespace ibmisc;
using namespace spsparse;
using namespace std::placeholders;





namespace icebin {
namespace modele {



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
    blitz::Array<double,1> const &foceanAOm,    // sparse indexing, 0-based
    blitz::Array<double,1> const &foceanAOp,    // sparse indexing, 0-based
    char invert = '+')
{
    SparseSetT &dimAOp(*ret.dims(1).sparse_set);

    // By looping over only things in IceBin ice sheet, we implicitly only
    // look at ice on ocean grid cells where both ModelE and IceBin have something.

    // Look only at ocean grid cells where BOTH ModelE and IcBin have ice
    for (int iAOp_d=0; iAOp_d < dimAOp.dense_extent(); ++iAOp_d) {
        long const iAO_s = dimAOp.to_sparse(iAOp_d);

        double const fcont_p = 1.0 - foceanAOp(iAO_s);
        double const fcont_m = 1.0 - foceanAOm(iAO_s);    // 0 or 1

        if (fcont_m == 0.0) continue;
        if (fcont_m != 1.0) (*icebin_error)(-1,
            "fcont_m must be 0 or 1");

        if (fcont_p == 0.0) continue;    // Can't extend ice that's not there

        // In Om space, cell only accounts for the fraction covered by continent.
        // When divided by wM (scaling), this will have the effect of magnifying
        // cells only partially covered by ocean (in PISM)
        ret.add({iAO_s, iAO_s}, invert == '-' ? fcont_p : 1./fcont_p);
    }
}
// ----------------------------------------------------------------
/** Helper function: clip a cell according to its index */
static bool dim_clip(SparseSetT const *dim, long index)
    { return dim->in_sparse(index); }


/** Matrix converts AO <- AA, based on Hntr.
Produces an UNSCALED (raw) matrix.

NOTES:
 1. Meant to be used with class MakeDenseEigenT.

@param hntrs
   Hntr grid definitions for {AO, AA}
@param wAO_d
   Weight (norm) of grid cells in AO. Dense indexing.
*/
static void raw_AOvAA(
    MakeDenseEigenT::AccumT &&ret,        // {dimA, dimO}
    std::array<HntrGrid const *, 2> hntrs,    // hntrAO, hntrAA
    blitz::Array<double,1> const &wAO_d)    // Size of each grid cell in AO
{
    auto &dimAO(*ret.dim(0).sparse_set);
    Hntr hntr_AOvAA(hntrs, 0);    // AOvAA
    hntr_AOvAA.matrix(std::move(ret),
        std::bind(&dim_clip, &dimAO, _1),
        *wAO_d);
}
// ----------------------------------------------------------------
/** Helper class for raw_EOvEA().
@see raw_EOvEA */
class RawEOvEA {
public:
    MakeDenseEigenT::AccumT ret;    // Final place for EOvEA
    GCMRegridder_ModelE const *gcmA;
    blitz::Array<double,1> const &wEO_d;

    RawEOvEA(
        MakeDenseEigenT::AccumT &&_ret,
        GCMRegridder_ModelE const *_gcmA,
        blitz::Array<double,1> const &_wEO_d)
    : ret(std::move(_ret)), gcmA(_gcmA), wEO_d(*_wEO_d) {}

    /** Called by Hntr::matrix() (AvO) */
    void add(std::array<long,2> index, double value)
    {
        // Ignore stray overlaps
        if (std::abs(value) < 1e-8) return;

        SparseSetT &dimEO(*ret.dim(0).sparse_set);
        long lAO_s = index[0];
        long lAA_s = index[1];

        // Iterate through all possible elevation classes for this gridcell pair
        for (int ihc=0; ihc<gcmA->nhc(); ++ihc) {
            long const lEO_s = gcmA->gcmO->indexingHC.tuple_to_index(
                std::array<long,2>{lAO_s,ihc});

            // Obtain wEO, size of the elevation grid cell
            if (!dimEO.in_dense(lEO_s)) continue;    // wEO==0 here
            int const lEO_d = dimEO.to_dense(lEO_s);
            double const weightEO = wEO_d(lEO_d);

            if (weightEO != 0) {
                int const lEA_s = gcmA->indexingHC.tuple_to_index(
                    std::array<long,2>{lAA_s,ihc});
                ret.add({lEO_s,lEA_s}, weightEO);
            }
        }
    }
};

/** Computes EOvEA (raw, unscaled).

Strategy: Group and sum EO elevation classes to produce EA elevation
   classes.  This is done by calling Hntr to discover which grid cells
   in AO correspond to each grid cell in AA.  Weights provided by Hntr
   are ignored, summing up wEO instead.

NOTES:

   1. This strategy assumes that every cell in AO is fully contained
      in a cell of AO.

   2. We should have:
               sum(EOvEA, 1) == wEA

      This can and should be checked; and if it does not work out, we
      should correct by multiplying by:
          EOvEA_corrected = EOvEA_raw * diag(sum(EOvEA,1,'-') .* wEA)
*/
static void raw_EOvEA(
    MakeDenseEigenT::AccumT &&ret,        // {dimEA, dimEO}; dimEO should not change here.
    std::array<HntrGrid const *, 2> hntrs,    // {hntrAO, hntrAA}
    SparseSetT const *dimEO,            // Used to clip in Hntr::matrix()
    GCMRegridder_ModelE const *gcmA,
    blitz::Array<double,1> *wEO_d)            // == EOvI.wM.  Dense indexing.
{
    // Call Hntr to generate AOvAA; and use that (above) to produce EOvEA
    Hntr hntr_AOvAA(hntrs, 0);    // dimB=A,  dimA=O
    RawEOvEA raw(std::move(ret), gcmA, wEO);
    hntr_AOvAA.matrix(raw,
        std::bind(&dim_clip, dimEO, _1));
}
// ----------------------------------------------------------------
// ========================================================================

/** Top-level subroutine produces the matrix EAmvIp or AAmvIp.

@param X
    Controls which matrix is produced:
      'E' -> EAmvIp
      'A' -> AAmvIp
@parm rmO
    Regrids between (AOp, EOp, Ip)
    (NOTE: rmO expects to see these matrices named "A", "E" and "I")
*/
static std::unique_ptr<WeightedSparse> compute_XAmvIp(
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &paramsA,
    char X,    // Either 'A' or 'E', determines the destination matrix.
    RegridMatrices const &rmO)
{
    TmpAlloc tmp;

    SparseSetT &dimXAm(*dims[0]);
    SparseSetT &dimIp(*dims[1]);
    SparseSetT dimEOm, dimEOp, dimOm;

    std::unique_ptr<WeightedSparse> ret(new WeightedSparse(dims, false));    // not conservative

    // ------------ Params for generating sub-matrices
    RegridMatrices::Params paramsO(paramsA);
    paramsO.scale = false;

    // ---------------- Take along the non-conservative weights

    // ----------- Determine universe for AOp / EOp
    // (and we will use the main matrix later)

    // We need AOpvIp for wAOP; used below.
    SparseSetT dimAOp;
    std::unique_ptr<WeightedSparse> AOpvIp(rmO.matrix("AvI", {&dimAOp, &dimIp}, paramsO));
    blitz::Array<double,1> const &wAOp(AOpvIp.mW);

    // dimAOm is a subset of dimAOp; we will use dimAOp to be sure.
    SparseSetT &dimAOm(dimAOp);

    std::unique_ptr<ConstUniverse> const_universe;

    // ----------- Compute weights for GCM ice sheet (wXOm)
    // wEOm = sEOmvOm * EOmvAOm * wc_AOmvAOp * wAOp
    //    or   wAOm = AOmvAOp * wAOp
    // wAOp

    // OmvOp
    EigenSparseMatrixT sc_AOmvAOp(MakeDenseEigenT(
        std::bind(&scaled_AOmvAOp, _1, foceanAOm, foceanAOp, '-'),
        {SparsifyTransform::TO_DENSE},
        {&dimAOm, &dimAOp}, '.').to_eigen());

    EigenColVectorT wAOm_e(scaled_AOmvAOp * map_eigen_colvector(wAOp));


    EigenColVectorT *wXOm_e;
    if (X == 'E') {
        // Get the universe for EOp / EOm
        SparseSetT dimEOp;
        const_universe.reset(new ConstUniverse({"dimIp"}, {&dimIp});
        std::unique_ptr<WeightedSparse> EOpvIp(rmO.matrix("EvI", {&dimEOp, &dimIp}, paramsO));
        const_universe.reset();        // Check that dimIp hasn't changed

        // dimEOm is a subset of dimEOp
        SparseSetT &dimEOm(dimEOp);

        // EOmvAOm: Repeat A values in E
        const_universe.reset(new ConstUniverse({"dimEOm", "dimAOm"}, {&dimEOm, &dimAOm});
        WeightedSparse EOmvAOm(rmO.matrix("EvA", {&dimEOm, &dimAOm}, paramsO));
        const_universe.reset();        // Check that dims didn't change
        auto sEOmvOm(1. / EOmvAOm->wM);

        // wEOm_e
        tmp.reset_ptr(wXOm_e, map_eigen_diagonal(sEOmvAOm) * *EOmvOm->M * wAOm_e);
    } else {
        wXOm_e = &wAOm_e;
    }

    // ---------------- Compute the main matrix
    // ---------------- XAmvIp = XAmvXOm * XOpvIp
    std::unique_ptr<EigenSparseMatrixT> XAmvXOm;    // Unscaled matrix
    if (X == 'E') {
        // Rename variables
        SparseSetT &dimEAm(dimXAm);
        blitz::Array<double,1> wEOm(to_blitz(wXOm_e));

        reset_ptr(XAmvXOm, MakeDenseEigenT(
            std::bind(&raw_EOvEA, _1, {&hntrAO, &hntrAA}, &dimEOm, to_blitz(*wXOm_e)),
            {SparsifyTransform::TO_DENSE_IGNORE_MISSING, SparsifyTransform::ADD_DENSE},
            {&dimEOm, &dimEAm}, 'T').to_eigen());
    } else {
        SparseSetT &dimAAm(dimXAm);

        // Actually AAmvAOm
        reset_ptr(XAmvXOm, MakeDenseEigenT(
            std::bind(&raw_AOvAA, _1, {&hntrAO, &hntrAA}, to_blitz(*wXOm_e)),
            {SparsifyTransform::TO_DENSE_IGNORE_MISSING, SparsifyTransform::ADD_DENSE},
            {&dimAOm, &dimAAm}, 'T').to_eigen());
    }

    // ---------- wEAm = [scale] * XAmvXOm * wEOm
    blitz::Array<double,1> sXAmvXOm(sum(XAmvXOm, 0, '-'));
    auto &wXAm_e(ret->tmp.make<EigenColVectorT>(
        map_eigen_diagonal(sXAmvXOm) * *XAmvXOm * *wXOm_e));

    // ----------- Put it all together (XAmvIp)
    ret->wM.reference(to_blitz(wXAm_e));
    if (paramsA.scale) {
        blitz::Array<double,1> sXAm(1. / to_blitz(wXAm_e));
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(sXAm) * *XAmvXOm * *XOpvIp->M));
    } else {
        ret->M.reset(new EigenSparseMatrixT(
            *XAmvXOm * *XOpvIp->M));
    }
    sum->Mw.reference(EOpvIp->Mw);
    return ret;
}
// ========================================================================
/** Creates Atmosphere grid from an existing Hntr-type Ocean grid. */
static std::unique_ptr<Grid> make_gridAA(Grid const *_gridAO)
{
    // -------- Check types on gridAO
    Grid_LonLat const *gridAO(dynamc_cast<Grid_LonLat *>(_gridAO));
    if (!gridAO) (*icebin_error)(-1,
        "make_gridAA() requires type Grid_LonLat");

    if (gridAO->north_pole != gridAO->south_pole) (*icebin_error)(-1,
        "north_pole=%d and south_pole=%d must match in gridAO",
        gridAO->north_pole, gridAO->south_pole);

    HntrGrid const *_hntrAO(&*gridAO->hntr);
    if (!_hntrAO) (*icebin_error)(-1,
        "make_gridAA() requires gridAO have a Hntr source");
    HntrGrid const &hntrAO(*_hntrAO);

    if ((hntrAO.im % 2 != 0) || (hntrAO.jm % 2 != 0)) (*icebin_error)(-1,
        "Ocean grid must have even number of gridcells for im and jm (vs. %d %d)",
        hntrAO.im, hntrAO.jm);

    // --------------------------
    // Define Atmosphere grid to be exactly twice the Ocean grid
    HntrGrid const hntrAA(hntrAO.im/2, hntrAO.jm/2, hntrAO.offi*0.5, hntrAO.dlat*2.);

    // Use hntr to figure out which grid cells should be realized in A, based
    // on realized grid cells in O
    SparseSetT dimA;
    Hntr hntrAOvA({&hntrAO, &hntrAA}, 0);
    hntrAOvA.matrix(
        accum::SparseSetAccum<SparseSetT,double,2>({nullptr, &dimA}),
        std::bind(&dim_clip, &gridAO->dim, _1));

    // ------- Produce a full gridAA, based on hntrAA and realized cells in dimA
    GridSpec_Hntr spec(hntrAA);
    spec.name = gridAO->name + "_A";
    spec.pole_caps = gridAO->north_pole;

    // Keep cells listed in dimA
    spec.spherical_clip = std::bind(&dim_clip, &dimA, _1);
    spec.points_in_side = 1;    // Keep it simple, we don't need this anyway
    spec.eq_rad = gridAO->eq_rad;

    std::unique_ptr<Grid_LonLat *> gridAA(new Grid_LonLat);
    spec.make_grid(*gridAA);
}

// ============================================================
// Methods for GCMRegridder_ModelE

GCMRegridder_ModelE::GCMRegridder_ModelE(
        std::unique_ptr<icebin::GCMRegridder> &&_gcmO,
        blitz::Array<double,1> const &_foceanAOm,
        blitz::Array<double,1> const &_foceanAOp) :
    : foceanAOm(reshape1(_foceanAOm)), foceanAOp(reshape1(foceanAOp)),
        gcmO(std::move(_gcmO))
{
    // Initialize baseclass members
    gridAA = make_gridAA(gcmO->gridAA);
    correctA = gcmO->correctA;
    hcdefs = gcmO->hcdefs;
    indexingHC = Indexing(
        {"O", "HC"},
        {0L, 0L},
        {gridAA->size(), gcmO->indexingHC[1].extent},
        {gcmO->indexingHC.indices()[0], gcmO->indexingHC.indices()[1]});
    indexingE = derive_indexingE(gridAA->indexing, indexingHC);
}

IceRegridder *GCMRegridder_ModelE::ice_regridder(std::string const &name) const
    { return gcmO->ice_regridder(name); }


RegridMatrices const GCMRegridder_ModelE::regrid_matrices(std::string const &ice_sheet_name) const
{
    IceModelRegridder const *ice_regridderO = gcmO->ice_regridder(ice_sheet_name);

    // Navigate the data structures
    Grid_LonLat const *gridAO = dynamic_cast<Grid_LonLat const *>(&*gcmO->gridAA);
    HntrGrid const &hntrAO(*gridAO->hntr);

    // Construct set of cells in A
    SparseSetT dimA;
    for (auto cell = gridAA->cells.begin(); cell != gridAA->cells.end(); ++cell)
        dimA.add_dense(cell->index);

    RegridMatrices rm;
    RegridMatrices const &rmO(rm.tmp.make<RegridMatrices>(
        ice_regridderO->regrid_matrices()));
    rm.add_regrid("EAmvIp", std::bind(&compute_EAmvIp, _1, _2,
        &dimA, &foceanAOm, &foceanAOp, &rmO));
    rm.add_regrid("AmvIp", std::bind(&compute_AmvIp, _1, _2,
        &dimA, &foceanAOm, &foceanAOp, &rmO));

    return rm;
}

}}    // namespace
