#include <spsparse/eigen.hpp>
#include <icebin/modele/GCMRegridder_ModelE.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/GridSpec_Hntr.hpp>

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
    MakeDenseEigenT::AccumT &ret,        // {dimAOm, dimAOp}
    blitz::Array<double,1> const &foceanAOp,    // sparse indexing, 0-based
    blitz::Array<double,1> const &foceanAOm,    // sparse indexing, 0-based
    char invert = '+')
{
    SparseSetT &dimAOp(*ret.dim(1).sparse_set);

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

    RawEOvEA(
        MakeDenseEigenT::AccumT &&_ret,
        GCMRegridder_ModelE const *_gcmA,
        blitz::Array<double,1> const &_wEO_d)
    : ret(std::move(_ret)), gcmA(_gcmA), wEO_d(_wEO_d) {}

    /** Called by Hntr::matrix() (AvO) */
    void add(std::array<long,2> index, double value)
    {
        // Ignore stray overlaps
        if (std::abs(value) < 1e-8) (*icebin_error)(-1,
            "Found a stray overlap; what should we do about it?");

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

   1. This strategy assumes that every cell in EO is fully contained
      in a cell of AO.

   2. We should have:
               sum(EOvEA, 1) == wEA

      This can and should be checked; and if it does not work out, we
      should correct by multiplying by:
          EOvEA_corrected = EOvEA_raw * diag(sum(EOvEA,1,'-') .* wEA)
*/
static void raw_EOvEA(
    MakeDenseEigenT::AccumT &ret,        // {dimEA, dimEO}; dimEO should not change here.
    std::array<HntrGrid const *, 2> hntrs,    // {hntrO, hntrA}
    SparseSetT const *dimAO,            // Used to clip in Hntr::matrix()
    GCMRegridder_ModelE const *gcmA,
    blitz::Array<double,1> &wEO_d)            // == EOvI.wM.  Dense indexing.
{
    // Call Hntr to generate AOvAA; and use that (above) to produce EOvEA
    Hntr hntr_AOvAA(hntrs, 0);    // dimB=A,  dimA=O
    RawEOvEA raw(std::move(ret), gcmA, wEO_d);

    hntr_AOvAA.overlap(raw, 1.0, DimClip(dimAO));
}


// ----------------------------------------------------------------
// ========================================================================
class ConstUniverse {
    std::vector<std::string> names;
    std::vector<SparseSetT *> dims;
    std::vector<int> extents;

public:
    ConstUniverse(
        std::vector<std::string> &&_names,
        std::vector<SparseSetT *> &&_dims) :
        names(std::move(_names)), dims(std::move(_dims))
    {
        if (names.size() != dims.size()) (*icebin_error)(-1,
            "names.size() and dims.size() must match");

        extents.reserve(dims.size());
        for (size_t i=0; i<dims.size(); ++i)
            extents.push_back(dims[i]->dense_extent());
    }

    ~ConstUniverse()
    {
        bool err = false;
        for (size_t i=0; i<dims.size(); ++i) {
            if (extents[i] != dims[i]->dense_extent()) {
                fprintf(stderr, "Dimension %s changed from %d to %d\n",
                    names[i].c_str(), extents[i], dims[i]->dense_extent());
                err = true;
            }
        }
        if (err) (*icebin_error)(-1,
            "At least one dimension changed");
    }
};


static Grid_LonLat const *cast_Grid_LonLat(Grid const *_gridO)
{
    // -------- Check types on gridO
    Grid_LonLat const *gridO(dynamic_cast<Grid_LonLat const *>(_gridO));
    if (!gridO) (*icebin_error)(-1,
        "make_gridA() requires type Grid_LonLat");

    if (gridO->north_pole != gridO->south_pole) (*icebin_error)(-1,
        "north_pole=%d and south_pole=%d must match in gridO",
        gridO->north_pole, gridO->south_pole);

    HntrGrid const *_hntrO(&*gridO->hntr);
    if (!_hntrO) (*icebin_error)(-1,
        "make_gridA() requires gridO have a Hntr source");

    return gridO;
}


static HntrGrid const make_hntrA(Grid_LonLat const *gridO)
{
    HntrGrid const &hntrO(*gridO->hntr);

    if ((hntrO.im % 2 != 0) || (hntrO.jm % 2 != 0)) (*icebin_error)(-1,
        "Ocean grid must have even number of gridcells for im and jm (vs. %d %d)",
        hntrO.im, hntrO.jm);

    // --------------------------
    // Define Atmosphere grid to be exactly twice the Ocean grid
    HntrGrid const hntrA(hntrO.im/2, hntrO.jm/2, hntrO.offi*0.5, hntrO.dlat*2.);

    return hntrA;
}


/** Creates Atmosphere grid from an existing Hntr-type Ocean grid. */
static std::unique_ptr<Grid> make_gridA(Grid_LonLat const *gridO)
{
printf("BEGIN make_gridA(%p)\n", gridO);

    HntrGrid const &hntrO(*gridO->hntr);
    HntrGrid const hntrA(make_hntrA(gridO));

    // Determine set of used grid cells in gridO
    SparseSetT dimO;
    for (auto cell=gridO->cells.begin(); cell != gridO->cells.end(); ++cell) {
        dimO.add_dense(cell->index);
    }

    // Use hntr to figure out which grid cells should be realized in A, based
    // on realized grid cells in O
    SparseSetT dimA;
    Hntr hntrOvA({&hntrO, &hntrA}, 0);

    // blitz::Array<double,1> wtO(dim_clip(dimO));

    accum::SparseSetAccum<SparseSetT,double,2> acc({nullptr, &dimA});
    hntrOvA.overlap(acc, 1.0, DimClip(&dimO));

    // ------- Produce a full gridA, based on hntrA and realized cells in dimA
    GridSpec_Hntr spec(hntrA);
    spec.name = gridO->name + "_A";
    spec.pole_caps = gridO->north_pole;

    // Keep cells listed in dimA
    spec.spherical_clip = std::bind(&dim_clip_fn, &dimA, _1);
    spec.points_in_side = 1;    // Keep it simple, we don't need this anyway
    spec.eq_rad = gridO->eq_rad;

    std::unique_ptr<Grid_LonLat> gridA(new Grid_LonLat);
    spec.make_grid(*gridA);

    std::unique_ptr<Grid> ret(gridA.release());

printf("END make_gridA(%p -> %p)\n", gridO, ret.get());
    return ret;
}

// ========================================================================
// --------------------------------------------------------
/** This exists mostly so we can test it from Python */
static std::unique_ptr<WeightedSparse> compute_AOmvAAm(
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &paramsA,
    GCMRegridder_ModelE const *gcmA,
    char transpose)
{
    auto &dimAOm(*dims[0]);
    auto &dimAAm(*dims[1]);


    // Obtain hntrA and hntrO for process below.
    Grid_LonLat const *_gcmO = cast_Grid_LonLat(&*gcmA->gcmO->gridA);
    Grid_LonLat const *_gcmA = cast_Grid_LonLat(&*gcmA->gridA);

    Hntr hntr_XOmvXAm(std::array<HntrGrid const *,2>{&*_gcmO->hntr, &*_gcmA->hntr});

    std::unique_ptr<WeightedSparse> ret(new WeightedSparse(dims, true));    // conservative
    reset_ptr(ret->M, MakeDenseEigenT(
        std::bind(&Hntr::overlap<MakeDenseEigenT::AccumT,DimClip>,
            &hntr_XOmvXAm, _1, _gcmA->eq_rad, DimClip(&dimAOm)),
        {SparsifyTransform::TO_DENSE_IGNORE_MISSING, SparsifyTransform::ADD_DENSE},
        {&dimAOm, &dimAAm}, transpose).to_eigen());

    // Ignore the scale parameter

    ret->wM.reference(sum(*ret->M, 0, '+'));
    ret->Mw.reference(sum(*ret->M, 1, '+'));

    return ret;
}

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
    GCMRegridder_ModelE const *gcmA,
    char X,    // Either 'A' or 'E', determines the destination matrix.
    double const eq_rad,    // Radius of the earth
    RegridMatrices const *rmO)
{
    TmpAlloc tmp;

    SparseSetT &dimXAm(*dims[0]);
    SparseSetT &dimIp(*dims[1]);
    SparseSetT dimEOm, dimEOp, dimOm;

    dimXAm.set_sparse_extent(X == 'A' ? gcmA->nA() : gcmA->nE());
    dimIp.set_sparse_extent(rmO->ice_regridder->nI());

    std::unique_ptr<WeightedSparse> ret(new WeightedSparse(dims, false));    // not conservative

    // ------------ Params for generating sub-matrices
    RegridMatrices::Params paramsO(paramsA);
    paramsO.scale = false;
    paramsO.correctA = true;

    // ---------------- Take along the non-conservative weights

    // ----------- Determine universe for AOp / EOp
    // (and we will use the main matrix later)

    // We need AOpvIp for wAOP; used below.
    SparseSetT dimAOp;
    std::unique_ptr<WeightedSparse> AOpvIp(rmO->matrix("AvI", {&dimAOp, &dimIp}, paramsO));
    blitz::Array<double,1> const &wAOp(AOpvIp->wM);

    // dimAOm is a subset of dimAOp; we will use dimAOp to be sure.
    SparseSetT &dimAOm(dimAOp);

    std::unique_ptr<ConstUniverse> const_universe;

    // ----------- Compute weights for GCM ice sheet (wXOm)
    // wEOm = sEOmvOm * EOmvAOm * wc_AOmvAOp * wAOp
    //    or   wAOm = AOmvAOp * wAOp
    // wAOp

    // OmvOp
    EigenSparseMatrixT sc_AOmvAOp(MakeDenseEigenT(
        std::bind(&scaled_AOmvAOp, _1, gcmA->foceanAOp, gcmA->foceanAOm, '+'),
        {SparsifyTransform::TO_DENSE},
        {&dimAOm, &dimAOp}, '.').to_eigen());

    EigenColVectorT wAOm_e(sc_AOmvAOp * map_eigen_colvector(wAOp));


    std::unique_ptr<EigenSparseMatrixT> XAmvXOm;    // Unscaled matrix

    // Obtain hntrA and hntrO for process below.
    HntrGrid const &hntrA(*cast_Grid_LonLat(&*gcmA->gridA)->hntr);
    HntrGrid const &hntrO(*cast_Grid_LonLat(&*gcmA->gcmO->gridA)->hntr);

    EigenColVectorT *wXOm_e;
    WeightedSparse *XOpvIp;
    if (X == 'E') {
        // Get the universe for EOp / EOm
        SparseSetT dimEOp;
        const_universe.reset(new ConstUniverse({"dimIp"}, {&dimIp}));
        WeightedSparse &EOpvIp(tmp.take(rmO->matrix("EvI", {&dimEOp, &dimIp}, paramsO)));

//        std::unique_ptr<WeightedSparse> EOpvIp();
        XOpvIp = &EOpvIp;

        const_universe.reset();        // Check that dimIp hasn't changed

        // dimEOm is a subset of dimEOp
        SparseSetT &dimEOm(dimEOp);

        // EOmvAOm: Repeat A values in E
        const_universe.reset(new ConstUniverse({"dimEOm", "dimAOm"}, {&dimEOm, &dimAOm}));
        std::unique_ptr<WeightedSparse> EOmvAOm(
            rmO->matrix("EvA", {&dimEOm, &dimAOm}, paramsO));
        const_universe.reset();        // Check that dims didn't change
        blitz::Array<double,1> sEOmvAOm(1. / EOmvAOm->wM);


        // wEOm_e
        tmp.take<EigenColVectorT>(wXOm_e,
            EigenColVectorT(map_eigen_diagonal(sEOmvAOm) * *EOmvAOm->M * wAOm_e));

        // ---------------- Compute the main matrix
        // ---------------- XAmvIp = XAmvXOm * XOpvIp
        // Rename variables
        SparseSetT &dimEAm(dimXAm);
        blitz::Array<double,1> wEOm(to_blitz(*wXOm_e));

        blitz::Array<double,1> wXOm(to_blitz(*wXOm_e));
        reset_ptr(XAmvXOm, MakeDenseEigenT(
            std::bind(&raw_EOvEA, _1,
                std::array<HntrGrid const *,2>{&hntrO, &hntrA},
                &dimEOm, gcmA, wXOm),
            {SparsifyTransform::TO_DENSE_IGNORE_MISSING, SparsifyTransform::ADD_DENSE},
            {&dimEOm, &dimEAm}, 'T').to_eigen());

    } else {    // X == 'A'
        XOpvIp = &*AOpvIp;
        wXOm_e = &wAOm_e;

        // ---------------- Compute the main matrix
        // ---------------- XAmvIp = XAmvXOm * XOpvIp

        SparseSetT &dimAAm(dimXAm);

        // Actually AOmvAAm
        Hntr hntr_XOmvXAm(std::array<HntrGrid const *,2>{&hntrO, &hntrA});
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
    blitz::Array<double,1> XAmvXOms(sum(*XAmvXOm, 1, '-'));    // Not unusual way we weight/scale here

    // The matrix *XAmvXOm * map_eigen_diagonal(XAmvXOms) should consist only of 1's and 0's
    auto &wXAm_e(ret->tmp.make<EigenColVectorT>(
        *XAmvXOm * map_eigen_diagonal(XAmvXOms) * *wXOm_e));

    // ----------- Put it all together (XAmvIp)
    ret->wM.reference(to_blitz(wXAm_e));
    if (paramsA.scale) {
        auto wXAm(to_blitz(wXAm_e));
        blitz::Array<double,1> sXAm(1. / wXAm);
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(sXAm) * *XAmvXOm * *XOpvIp->M));
    } else {
        ret->M.reset(new EigenSparseMatrixT(
            *XAmvXOm * *XOpvIp->M));
    }
    ret->Mw.reference(XOpvIp->Mw);
    return ret;
}
// ========================================================================

// ============================================================
// Methods for GCMRegridder_ModelE

GCMRegridder_ModelE::GCMRegridder_ModelE(
        std::shared_ptr<icebin::GCMRegridder> const &_gcmO)
    : gcmO(_gcmO)
{
    // Initialize baseclass members
    gridA = make_gridA(cast_Grid_LonLat(&*gcmO->gridA));
    correctA = gcmO->correctA;
    hcdefs = gcmO->hcdefs;
    indexingHC = Indexing(
        {"O", "HC"},
        {0L, 0L},
        {gridA->ndata(), gcmO->indexingHC[1].extent},
        {gcmO->indexingHC.indices()[0], gcmO->indexingHC.indices()[1]});
    indexingE = derive_indexingE(gridA->indexing, indexingHC);
}

void GCMRegridder_ModelE::set_focean(
        blitz::Array<double,1> &_foceanAOp,
        blitz::Array<double,1> &_foceanAOm)
{
    // Take shared referenes, to prevent de-allocation
    this->_foceanAOp.reference(_foceanAOp);
    this->_foceanAOm.reference(_foceanAOm);

    // Use reshaped versions, which don't share a reference
    foceanAOp.reference(reshape1(this->_foceanAOp));
    foceanAOm.reference(reshape1(this->_foceanAOm));    // set lbound=0
}

IceRegridder *GCMRegridder_ModelE::ice_regridder(std::string const &name) const
    { return gcmO->ice_regridder(name); }


RegridMatrices GCMRegridder_ModelE::regrid_matrices(std::string const &ice_sheet_name) const
{
    Grid_LonLat const *gridO(cast_Grid_LonLat(&*gcmO->gridA));

    RegridMatrices _rmO(gcmO->regrid_matrices(ice_sheet_name));
    RegridMatrices rm(_rmO.ice_regridder);
    RegridMatrices &rmO(rm.tmp.take<RegridMatrices>(std::move(_rmO)));

    rm.add_regrid("EAmvIp", std::bind(&compute_XAmvIp, _1, _2,
        this, 'E', gridO->eq_rad, &rmO));
    rm.add_regrid("AAmvIp", std::bind(&compute_XAmvIp, _1, _2,
        this, 'A', gridO->eq_rad, &rmO));

    // For testing
    rm.add_regrid("AOmvAAm", std::bind(&compute_AOmvAAm, _1, _2,
        this, '.'));
    rm.add_regrid("AAmvAOm", std::bind(&compute_AOmvAAm, _1, _2,
        this, 'T'));

    return std::move(rm);
}

}}    // namespace
