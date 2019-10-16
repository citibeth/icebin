#include <icebin/modele/topo.hpp>
#include <spsparse/eigen.hpp>
#include <ibmisc/stdio.hpp>
#include <ibmisc/linear/compressed.hpp>
#include <icebin/modele/GCMRegridder_ModelE.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/merge_topo.hpp>
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
// ----------------------------------------------------------------



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
    GridSpec_LonLat const &specO(gcmA->specO());
    GridSpec_LonLat const &specA(gcmA->specA());

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
// ========================================================================
// ------------------------------------------------------
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
        blitz::Array<double,1> const &foceanAOp,
        blitz::Array<double,1> const &foceanAOm,
        char X,    // Either 'A' or 'E', determines the destination matrix.
        double const eq_rad,    // Radius of the earth
        RegridMatrices_Dynamic const *rmO);
};    // class ComputeXAmvIp_Helper

ComputeXAmvIp_Helper::ComputeXAmvIp_Helper(
    std::array<SparseSetT *,2> dims,
    RegridParams const &paramsA,
    GCMRegridder_ModelE const *gcmA,
    blitz::Array<double,1> const &foceanAOp,
    blitz::Array<double,1> const &foceanAOm,
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

    EigenColVectorT wAOm_e(compute_wAOm(foceanAOp, foceanAOm, wAOp, dimAOp, dimAOm));

    HntrSpec const &hntrA(cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr);
    HntrSpec const &hntrO(cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr);

    dimXAm.set_sparse_extent(X == 'A' ? gcmA->nA() : gcmA->nE());
    dimIp.set_sparse_extent(rmO->ice_regridder->nI());

    // Ensure that operations don't change universes.  Sanity check.
    std::unique_ptr<ConstUniverseT> const_universe;

    if (X == 'E') {
        RegridParams paramsO(paramsA);
        paramsO.scale = false;
        paramsO.correctA = false;

        // Get the universe for EOp / EOm
        const_universe.reset(new ConstUniverseT({"dimIp"}, {&dimIp}));
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
    blitz::Array<double,1> const &foceanAOp,    // gcmA->foceanAOp
    blitz::Array<double,1> const &foceanAOm,    // gcmA->foceanAOm
    char X,    // Either 'A' or 'E', determines the destination matrix.
    double const eq_rad,    // Radius of the earth
    RegridMatrices_Dynamic const *rmO)
{
    std::string dimXA(strprintf("dim%cA", X));
    std::unique_ptr<linear::Weighted_Eigen> ret(new linear::Weighted_Eigen(dims, false));    // not conservative

    // Do the legwork, and fish out things from the results that we need
    ComputeXAmvIp_Helper hh(dims, paramsA, gcmA, foceanAOp, foceanAOm, X, eq_rad, rmO);
    ret->tmp = std::move(hh.tmp);
    auto &XOpvIp(hh.XOpvIp);
    auto &wXAm_e(*hh.wXAm_e);
    auto &XAmvXOm(hh.XAmvXOm);

    auto &dimXOm(X=='E' ? hh.dimEOm : hh.dimAOm);
    auto &dimXOp(X=='E' ? hh.dimEOp : hh.dimAOp);

//printf("XOpvIp: %ld\n", (long)(XOpvIp->M->nonZeros()));
//printf("XAmvXOm: %ld\n", (long)(XAmvXOm->nonZeros()));

    // ----------- Put it all together (XAmvIp)
    blitz::Array<double,1> sXOpvIp(1. / XOpvIp->wM);
    blitz::Array<double,1> wXAm(to_blitz(wXAm_e));
    ret->wM.reference(wXAm);    // wAAm_e or wEAm_e
    blitz::Array<double,1> sXAm(sum(*XAmvXOm, 0, '-'));
    if (paramsA.scale) {
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(sXAm) *  *XAmvXOm *
            crop_mvp(dimXOm, dimXOp, 0,
                map_eigen_diagonal(sXOpvIp) * *XOpvIp->M)));
    } else {
        // Apply fixes done for scaled version (2017-10-15) to unscaled
        // see 5e0f7d79.  Ensure that M_unscaled = wM * M_scaled
        blitz::Array<double,1> scale(sXAm.extent(0));
        for (int i=0; i<sXAm.extent(0); ++i) scale(i) = wXAm(i) * sXAm(i);
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(scale) *  *XAmvXOm *
            crop_mvp(dimXOm, dimXOp, 0,
                map_eigen_diagonal(sXOpvIp) * *XOpvIp->M)));
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
    blitz::Array<double,1> const &foceanAOp,    // gcmA->foceanAOp
    blitz::Array<double,1> const &foceanAOm,    // gcmA->foceanAOm
    char X,    // Either 'A' or 'E', determines the destination matrix.
    double const eq_rad,    // Radius of the earth
    RegridMatrices_Dynamic const *rmO)
{
    std::string dimXA(strprintf("dim%cA", X));
    std::unique_ptr<linear::Weighted_Eigen> ret(new linear::Weighted_Eigen(dims, false));    // not conservative

    // Do the legwork, and fish out things from the results that we need
    ComputeXAmvIp_Helper hh({dims[1], dims[0]}, paramsA, gcmA, foceanAOp, foceanAOm, X, eq_rad, rmO);
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
    std::shared_ptr<icebin::GCMRegridder_Standard> const &_gcmO)
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

    // Read number of global EC's out of global EC matrix file.
    if (global_ecO == "") {
        global_hcdefs.clear();
    } else {
printf("Opening global_ecO = %s\n", global_ecO.c_str());
        NcIO ncio(global_ecO, 'r');
        ncio_vector(ncio, global_hcdefs, true, "hcdefs", "double",
            get_or_add_dims(ncio, {"nhc"}, {_hcdefs.size()} ));
    }

    // Combine local and global hcdefs
    _hcdefs.clear();
    for (double elev : gcmO->hcdefs()) _hcdefs.push_back(elev);
    for (double elev : global_hcdefs) _hcdefs.push_back(elev);

    // Read EOpvAOp_base from global_ec file
    // Read metadata and global EOpvAOp matrix (from output of global_ec.cpp)
    if (global_ecO != "") {
        {NcIO ncio(global_ecO, 'r');
            // metaO.ncio(ncio);   // no metaO in this class
            EOpvAOp_base.ncio(ncio, "EvO.M");
        }
    }

}




std::unique_ptr<RegridMatrices_Dynamic> GCMRegridder_ModelE::regrid_matrices(
    int sheet_index,
    blitz::Array<double,1> const &foceanAOp,
    blitz::Array<double,1> const &foceanAOm,
    blitz::Array<double,1> const &elevmaskI,
    RegridParams const &params) const
{
    IceRegridder const *regridder = &*ice_regridders()[sheet_index];
    std::unique_ptr<RegridMatrices_Dynamic> rm(new RegridMatrices_Dynamic(regridder, params));
    GridSpec_LonLat const &specO(this->specO());

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
        this, foceanAOp, foceanAOm,
        'E', specO.eq_rad, &rmO));
    rm->add_regrid("AAmvIp", std::bind(&compute_XAmvIp, _1, _2,
        this, foceanAOp, foceanAOm,
        'A', specO.eq_rad, &rmO));

// We never need this.
//    rm->add_regrid("AAmvEAm", std::bind(&compute_AAmvEAm_rmO, _1, _2,
//        this, specO.eq_rad, &rmO));


    rm->add_regrid("IpvEAm", std::bind(&compute_IpvXAm, _1, _2,
        this, foceanAOp, foceanAOm,
        'E', specO.eq_rad, &rmO));
    rm->add_regrid("IpvAAm", std::bind(&compute_IpvXAm, _1, _2,
        this, foceanAOp, foceanAOm,
        'A', specO.eq_rad, &rmO));


    // --------------- Simply-named aliases
    rm->add_regrid("EvI", std::bind(&compute_XAmvIp, _1, _2,
        this, foceanAOp, foceanAOm,
        'E', specO.eq_rad, &rmO));
    rm->add_regrid("AvI", std::bind(&compute_XAmvIp, _1, _2,
        this, foceanAOp, foceanAOm,
        'A', specO.eq_rad, &rmO));


//    rm->add_regrid("AvE", std::bind(&compute_AAmvEAm_rmO, _1, _2,
//        this, specO.eq_rad, &rmO));


    rm->add_regrid("IvE", std::bind(&compute_IpvXAm, _1, _2,
        this, foceanAOp, foceanAOm,
        'E', specO.eq_rad, &rmO));
    rm->add_regrid("IvA", std::bind(&compute_IpvXAm, _1, _2,
        this, foceanAOp, foceanAOm,
        'A', specO.eq_rad, &rmO));

    // For testing
    rm->add_regrid("AOmvAAm", std::bind(&compute_AOmvAAm, _1, _2,
        this, '.'));
    rm->add_regrid("AAmvAOm", std::bind(&compute_AOmvAAm, _1, _2,
        this, 'T'));

    return rm;
}
// -----------------------------------------------------------------------
/** Computes global AvE, including any base ice, etc.
    @param emI_lands One emI_land array per ice sheet (elevation on continent, NaN in ocean).
    @param emI_ices One emI_ice array per ice sheet (elevation on ice, NaN off ice).
    @param params Parameters to use in generating regridding matrices.
        Should be RegridParams(true, true, {0,0,0}) to give conservative matrix.
    @param offsetE Offset (in sparse E space) added to base EC indices */
linear::Weighted_Tuple GCMRegridder_ModelE::global_AvE(
    std::vector<blitz::Array<double,1>> const &emI_lands,
    std::vector<blitz::Array<double,1>> const &emI_ices,
    blitz::Array<double,1> const &foceanAOp,
    blitz::Array<double,1> const &foceanAOm,
    bool scale,
    // ----------- Output vars
    long &offsetE) const
{
    GridSpec_LonLat const &specO(this->specO());

    // --------------------- Compute EOpvAOp (merged global + local ice)

    // struct EOpvAOpResult:
    //     SparseSetT dimEOp;    // dimEOp is set and returned; dimAOp is appended
    //     std::unique_ptr<EigenSparseMatrixT> EOpvAOp;
    //     std::vector<double> hcdefs; // OUT:  Elev class definitions for merged ice
    //     ibmisc::Indexing indexingHC;
    //     std::vector<uint16_t> underice_hc;

    std::vector<std::string> errors;
    SparseSetT dimAOp;
    EOpvAOpResult eam(compute_EOpvAOp_merged(
        dimAOp, EOpvAOp_base,
        RegridParams(false, false, {0.,0.,0.}),  // (scale, correctA, sigma)
        &*gcmO, specO.eq_rad, emI_ices,
        true, true,    // use_global_ice=t, use_local_ice=t
        hcdefs(), indexingHC, false, errors));
    offsetE = eam.offsetE;   // Return offsetE value

    // Print sanity check errors to STDERR
    if (errors.size() > 0) {
        for (std::string const &err : errors) fprintf(stderr, "ERROR: %s\n", err.c_str());
        (*icebin_error)(-1, "GCMRegridder_ModelE: Problems detected creating merged global_ecO file, stopping now.");
    }


    // ----------------- Compute AAMvEAm
    auto wAOp(sum(*eam.EOpvAOp, 1, '+'));
    SparseSetT dimAAm, dimEAm;
    return _compute_AAmvEAm(
        scale,
        specO.eq_rad,
        cast_GridSpec_LonLat(*this->gcmO->agridA.spec).hntr,
        cast_GridSpec_LonLat(*this->agridA.spec).hntr,
        this->gcmO->indexingHC,
        this->indexingHC,
        foceanAOp,
        foceanAOm,
        *eam.EOpvAOp, eam.dimEOp, dimAOp, wAOp);


#if 0
    // ---------------- Convert to sparse indexing as linear::Weighted_Tuple
    linear::Weighted_Tuple AvE_g;

    spcopy(
        accum::to_sparse(make_array(&dimAAm),
        accum::ref(AvE_g.wM)),    // Output
        AAmvEAm.wM);   // Input

    spcopy(
        accum::to_sparse(make_array(E1vI.dims[0], &dimE0),
        accum::ref(AvE_g.M)),
        *AAmvEAm.M);

    spcopy(
        accum::to_sparse(make_array(&dimE0),
        accum::ref(AvE_g.Mw)),
        AAmvEAm.Mw);

    return AvE_g;
#endif

}

linear::Weighted_Tuple GCMRegridder_ModelE::global_unscaled_E1vE0(
    std::vector<linear::Weighted_Eigen *> const &E1vIs_unscaled, // State var set in IceCoupler::couple(); _nc = no correctA (RegridParam) UNSCALED matrix
    std::vector<EigenSparseMatrixT *> const &IvE0s, // State var set in IceCoupler::couple()  SCALED matrix
    std::vector<SparseSetT *> const &dimE0s) const    // dimE0 accompanying IvE0
{
    linear::Weighted_Tuple E1vE0_g;    // _g = global

    // ---------- Compute each E1vE0 and merge...
    for (size_t sheet_index=0; sheet_index < ice_regridders().size(); ++sheet_index) {
        auto &E1vI_unscaled(*E1vIs_unscaled[sheet_index]);
        EigenSparseMatrixT &IvE0(*IvE0s[sheet_index]);
        auto &dimE0(*dimE0s[sheet_index]);

        // Don't do this on the first round, since we don't yet have an IvE0
        EigenSparseMatrixT E1vE0(*E1vI_unscaled.M * IvE0);    // UNSCALED
        spcopy(
            accum::to_sparse(make_array(E1vI_unscaled.dims[0]),
            accum::ref(E1vE0_g.wM)),    // Output
            E1vI_unscaled.wM);   // Input

        spcopy(
            accum::to_sparse(make_array(E1vI_unscaled.dims[0], &dimE0),
            accum::ref(E1vE0_g.M)),
            E1vE0);

#if 0    // Not needed
        spcopy(
            accum::to_sparse(make_array(&dimE0),
            accum::ref(E1vE0_g.Mw)),
            IvE0.Mw);
#endif
    }    // Flush accumulators

    E1vE0_g.set_shape(std::array<long,2>{nE(), nE()});
    return E1vE0_g;
}

// ==============================================================

void GCMRegridder_WrapE::_init(std::unique_ptr<GCMRegridder_ModelE> &&_gcmA)
{
    gcmA = std::move(_gcmA);
    _ice_regridders = &gcmA->ice_regridders();

    // Copy assorted stuff over
    // agridA = gcmA->agridA;    // hopefully we don't have to copy this
    correctA = gcmA->correctA;
    indexingHC = gcmA->indexingHC;
    indexingE = gcmA->indexingE;
    _hcdefs = gcmA->_hcdefs;
}


GCMRegridder_WrapE::GCMRegridder_WrapE(std::unique_ptr<GCMRegridder_ModelE> &&_gcmA)
{
    _init(std::move(_gcmA));

    auto const nO = gcmA->gcmO->nA();
    foceanOp.reference(blitz::Array<double,1>(nO));
    foceanOm.reference(blitz::Array<double,1>(nO));
}

GCMRegridder_WrapE::GCMRegridder_WrapE(
    std::unique_ptr<GCMRegridder_ModelE> &&_gcmA,
    blitz::Array<double,1> _foceanOp,
    blitz::Array<double,1> _foceanOm)
: foceanOp(_foceanOp), foceanOm(foceanOm)
{
    _init(std::move(_gcmA));
}



}}    // namespace
