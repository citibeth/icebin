#include <spsparse/eigen.hpp>
#include <ibmisc/stdio.hpp>
#include <icebin/modele/GCMRegridder_ModelE.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/gridgen/GridGen_LonLat.hpp>

using namespace icebin;
using namespace ibmisc;
using namespace spsparse;
using namespace std::placeholders;

namespace spsparse {
namespace accum {

template<class IndexT, int RANK>
class Null {
{
public:
    static const int rank = RANK;
    typedef IndexT index_type;

    void set_shape(std::array<long, rank> const &_shape)
        { }

    template<class ValueT>
    void add(
        std::array<index_type, rank> const &index,
        ValueT const &val)
        {  }
};

}}    // namespace



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

int n=0;

    RawEOvEA(
        MakeDenseEigenT::AccumT &&_ret,
        GCMRegridder_ModelE const *_gcmA,
        blitz::Array<double,1> const &_wEO_d)
    : ret(std::move(_ret)), gcmA(_gcmA), wEO_d(_wEO_d) {}

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
        for (int ihc=0; ihc<gcmA->nhc(); ++ihc) {
            long const lEO_s = gcmA->gcmO->indexingHC.tuple_to_index(
                std::array<long,2>{lAO_s,ihc});

            // Obtain wEO, size of the elevation grid cell
            if (!dimEO.in_sparse(lEO_s)) continue;    // wEO==0 here
            int const lEO_d = dimEO.to_dense(lEO_s);
            double const weightEO = wEO_d(lEO_d);
//if (n < 15) printf("      ihc=%d, lEO_s = %ld   lEO_d=%d  weightEO=%g\n", ihc, lEO_s, lEO_d, weightEO);

            if (weightEO != 0) {
                int const lEA_s = gcmA->indexingHC.tuple_to_index(
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
    GCMRegridder_ModelE const *gcmA,
    blitz::Array<double,1> &wEO_d)            // == EOvI.wM.  Dense indexing.
{
    // Call Hntr to generate AOvAA; and use that (above) to produce EOvEA
    Hntr hntr_AOvAA(17.17, hntrO, hntrA, 0);    // dimB=A,  dimA=O

    hntr_AOvAA.overlap<RawEOvEA, DimClip>(
        RawEOvEA(std::move(ret), gcmA, wEO_d),
        eq_rad, DimClip(dimAO));
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


HntrSpec make_hntrA(HntrSpec const &hntrO)
{
    if ((hntrO.im % 2 != 0) || (hntrO.jm % 2 != 0)) (*icebin_error)(-1,
        "Ocean grid must have even number of gridcells for im and jm (vs. %d %d)",
        hntrO.im, hntrO.jm);

    // --------------------------
    // Define Atmosphere grid to be exactly twice the Ocean grid
    return HntrSpec(hntrO.im/2, hntrO.jm/2, hntrO.offi*0.5, hntrO.dlat*2.);
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

// ------------------------------------------------------
/** Computes the intermediate value wAOm, which is needed by more than one
    top-level regrid generator.  Computes some other values along with wAOm as well. */
class Compute_wAOm {
public:
    RegridParams paramsO;   // Parameters used with rmO
    SparseSetT dimAOp; // Includes only AO grid cells that interact with an ice sheet.
    SparseSetT dimAOm; // Includes only AO grid cells that interact with an ice sheet AND are not rounded to ocean
    EigenColVectorT wAOm_e;

    /** Parameters come from top-level regridding subroutine */
    Compute_wAOm(
        SparseSetT &dimIp,
        RegridParams const &paramsA,
        GCMRegridder_ModelE const *gcmA,
        RegridMatrices_Dynamic const *rmO);
};

Compute_wAOm::Compute_wAOm(
    SparseSetT &dimIp,
    RegridParams const &paramsA,
    GCMRegridder_ModelE const *gcmA,
    RegridMatrices_Dynamic const *rmO)
: paramsO(paramsA)
{
    // ------------ Params for generating sub-matrices
    paramsO.scale = false;
    paramsO.correctA = true;

    // We need AOpvIp for wAOP; used below.
printf("dimIp.dense_extent() = %d (sparse=%d)\n", dimIp.dense_extent(), dimIp.sparse_extent());
    std::unique_ptr<linear::Weighted_Eigen> AOpvIp_corrected(rmO->matrix_d(
        "AvI", {&dimAOp, &dimIp}, paramsO));
    blitz::Array<double,1> const &wAOp(AOpvIp_corrected->wM);
printf("dimIp.dense_extent() = %d (sparse=%d)\n", dimIp.dense_extent(), dimIp.sparse_extent());


printf("|AOpvIp_corrected|=%ld\n", (long)AOpvIp_corrected->M->nonZeros());
int n=0;
for (auto ii(begin(*AOpvIp_corrected->M)); ii != end(*AOpvIp_corrected->M); ++ii, ++n) {
    printf("   AOpvIp_corrected(%d %d) = %g\n", dimAOp.to_sparse(ii->row()), dimIp.to_sparse(ii->col()), ii->value());
    if (n > 15) break;
}

printf("|wAOp|1=%d\n", wAOp.extent(0));
for (int i=0; i<15; ++i) {
    printf("    wAOp(%d)=%g\n", i, wAOp(i));
}


    // ----------- Compute dimAOm properly
    // dimAOm is a subset of dimAOp.  Remove points in dimAOp that are ocean.
    for (int iAOp_d=0; iAOp_d < dimAOp.dense_extent(); ++iAOp_d) {
        auto const iAO_s(dimAOp.to_sparse(iAOp_d));
        if (gcmA->foceanAOm(iAO_s) == 0) dimAOm.add_dense(iAO_s);
    }
    printf("|dimAOm|=%ld\n", (long)dimAOm.dense_extent());
    printf("|dimAOp|=%ld\n", (long)dimAOp.dense_extent());
#if 0
    HntrSpec const &hspecA(cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr);
    HntrSpec const &hspecO(cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr);
    blitz::Array<double, 1> WTO(const_array(blitz::shape(hspecO.jm * hspecO.im), 1.0));
    Hntr hntr_AAmvAOm(17.17, hspecA, hspecO);
    blitz::Array<double,1> foceanAAm(hntr_AAmvAOm.regrid(WTO, foceanAOm));
#endif


    // OmvOp
    EigenSparseMatrixT sc_AOmvAOp(MakeDenseEigenT(
        std::bind(&scaled_AOmvAOp, _1, gcmA->foceanAOp, gcmA->foceanAOm, '+'),
        {SparsifyTransform::TO_DENSE_IGNORE_MISSING},
        {&dimAOm, &dimAOp}, '.').to_eigen());

    wAOm_e = sc_AOmvAOp * map_eigen_colvector(wAOp);

printf("|wAOm_e|=%ld   |wAOp|=%ld\n", wAOm_e.rows()*wAOm_e.cols(), wAOp.extent(0));
for (int i=0; i<15; ++i) {
    printf("    wAOp(%d)=%g   wAOm_e(%d)=%g\n", i, wAOp(i), i, wAOm_e(i));
}

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
static std::unique_ptr<linear::Weighted_Eigen> compute_AAmvEAm(
    std::array<SparseSetT *,2> dims,
    RegridParams const &paramsA,
    GCMRegridder_ModelE const *gcmA,
    double const eq_rad,    // Radius of the earth
    RegridMatrices_Dynamic const *rmO)
{
    SparseSetT &dimAAm(*dims[0]);
    SparseSetT &dimEAm(*dims[1]);
    SparseSetT dimIp;

    // Must set sparse_extent
    dimAAm.set_sparse_extent(gcmA->nA());
    dimEAm.set_sparse_extent(gcmA->nE());

    std::unique_ptr<linear::Weighted_Eigen> ret(new linear::Weighted_Eigen(dims, false));    // not conservative

    // ------------- Compute wAOm and related quantities
    Compute_wAOm c1(dimIp, paramsA, gcmA, rmO);
    auto &paramsO(c1.paramsO);
    auto &dimAOp(c1.dimAOp);
    auto &dimAOm(c1.dimAOm);
    auto &wAOm_e(c1.wAOm_e);

    // ------------- Compute wAAm and AAmvAOm
    // Obtain hntrA and hntrO for process below.
    HntrSpec const &hntrA(cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr);
    HntrSpec const &hntrO(cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr);

    if (!hntrA.is_set()) (*icebin_error)(-1, "hntrA must be set");
    if (!hntrO.is_set()) (*icebin_error)(-1, "hntrO must be set");

/*
AOmvAAm
The problem is.... on A gridcells where SOME of the O grid cells are used and SOME are not, weighting is not being taken into account properly.  This de-weights fhc (AvE) by a fraction within that gridcell.

DimClip here needs dimAOm.  Instead, it is being given dimAOp, which is a superset of dimAOm (see comment in Compute_wAOm).  The result is hntr is being scaled as if AOp, but with values for AOm.  Thus the under-scaling of A grid cells containing P grid cells that have been zeroed out.
*/
    Hntr hntr_AOmvAAm(17.17, hntrO, hntrA);
    EigenSparseMatrixT AAmvAOm(MakeDenseEigenT(
        std::bind(&Hntr::overlap<MakeDenseEigenT::AccumT,DimClip>,
            &hntr_AOmvAAm, _1, eq_rad, DimClip(&dimAOm)),
        {SparsifyTransform::TO_DENSE_IGNORE_MISSING, SparsifyTransform::ADD_DENSE},
        {&dimAOm, &dimAAm}, 'T').to_eigen());

    blitz::Array<double,1> AAmvAOms(sum(AAmvAOm, 1, '-'));    // Note unusual way we weight/scale here
    auto &wAAm_e(ret->tmp.make<EigenColVectorT>(
        AAmvAOm * map_eigen_diagonal(AAmvAOms) * wAOm_e));


    // ------------ Compute AOmvEOm
    SparseSetT dimEOm;
    std::unique_ptr<linear::Weighted_Eigen> AOmvEOm(    // unscaled
        rmO->matrix_d("AvE", {&dimAOm, &dimEOm}, paramsO));
    auto sAOmvEOm(sum(*AOmvEOm->M, 0, '-'));
===> Why not just compute EvA here?  (Because we also need AvE below...)

    auto EOmvAOm(AOmvEOm->M->transpose());
    auto &EOmvAOms(sAOmvEOm);
    EigenColVectorT wEOm_e(EOmvAOm * map_eigen_diagonal(EOmvAOms) * wAOm_e);
    auto wEOm(to_blitz(wEOm_e));


    // ------------ Compute EOmvEAm
    EigenSparseMatrixT EOmvEAm(MakeDenseEigenT(    // TODO: Call this EAvEO, since it's the same 'm' or 'p'
        std::bind(&raw_EOvEA, _1,
            hntrO, hntrA,
            eq_rad, &dimAOm, gcmA, wEOm),
        {SparsifyTransform::TO_DENSE_IGNORE_MISSING, SparsifyTransform::ADD_DENSE},
        {&dimEOm, &dimEAm}, '.').to_eigen());

    auto EAmvEOm(EOmvEAm.transpose());

    // ------------ Compute wEAm
    auto EAmvEOms(sum(EOmvEAm, 0, '-'));
    auto &sEOmvEAm(EAmvEOms);
    auto &wEAm_e(ret->tmp.make<EigenColVectorT>(
        EAmvEOm * map_eigen_diagonal(EAmvEOms) * wEOm_e));


    // ------------- Put it all together
    ret->wM.reference(to_blitz(wAAm_e));
    blitz::Array<double,1> sAAmvAOm(sum(AAmvAOm, 0, '-'));
    if (paramsA.scale) {
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(sAAmvAOm) * AAmvAOm *    // Works on whole cells, not just ice sheet; so we need to normalize by sAAmvAOm, not sAAm
            map_eigen_diagonal(sAOmvEOm) * *AOmvEOm->M *
            map_eigen_diagonal(sEOmvEAm) * EOmvEAm));
    } else {
        blitz::Array<double,1> scale(ret->wM * sAAmvAOm);
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(scale) *
            AAmvAOm *    // Works on whole cells, not just ice sheet; so we need to normalize by sAAmvAOm, not sAAm
            map_eigen_diagonal(sAOmvEOm) * *AOmvEOm->M *
            map_eigen_diagonal(sEOmvEAm) * EOmvEAm));
    }
    ret->Mw.reference(to_blitz(wEAm_e));

    return ret;
}
// ========================================================================
static std::unique_ptr<linear::Weighted_Eigen> compute_EOmvAOm_unscaled(
    SparseSetT *dimAOm_ptr,        // const; pre-computed, should not change
    SparseSetT &dimEOp,
    SparseSetT &dimAOp,
//    std::array<SparseSetT *,2> dims,  // dimEOm to be computed; dimAOm is given
    RegridMatrices_Dynamic const *rmO,
    RegridParams paramsO)   // copy
{
    paramsO.scaled = false;

    std::unique_ptr<linear::Weighted_Eigen> ret(
        new linear::Weighted_Eigen(true));    // conservative
    ret->dims[0] = ret->tmp.newptr<SparseSetT>();
    ret->dims[1] = dimAOm_ptr;

    // EOpvAOp: Repeat A values in E
    const_universe.reset(new ConstUniverse({"dimEOp", "dimAOp"}, {&dimEOp, &dimAOp}));
    std::unique_ptr<linear::Weighted_Eigen> EOpvAOp(
        rmO->matrix_d("EvA", {&dimEOp, &dimAOp}, paramsO));

    const_universe.reset();        // Check that dims didn't change

    // Convert to EOmvAOm by removing cells not in dimAOm
    blitz::Array<double,1> EOmvAOms(dimAOm.dense_extent());
    for (auto ii=EOpvAOp->M->begin(); ii != EOpvAOp->M->end(); ++ii) {
        auto const iAOp_d = ii->index(1);
        auto const iAOp_s = dimAOp.to_sparse(iAOp_d);
        int iAOm_d;
        if (dimAOm.to_dense_ignore_missing(iAOp_s, iAOm_d)) {
            auto const iEOp_d = ii->index(1);
            auto const iEO_s = dimEOp.to_sparse(iEOp_d);
            auto const iEOm_d = dimEOm.add_dense(iEO_s);  // add to dimEOp
            EOmvAOm_tl.add({iEOm_d, iAOm_d}, ii->value()); // add to EOmvAOm
        }
    }

    EOmvAOm_tl.set_shape(std::array<long,2>{});
    ret->M.reset(new EigenSparseMatrixT(dimEOm.dense_extent(), dimAOm.dense_extent()));
    ret->M->setFromTriplets(EOmvAOm_tl.begin(), EOmvAOm_tl.end());

    ret->wM.reference(sum(*ret->M, 0, '+'));
    ret->Mw.reference(sum(*ret->M, 1, '+'));

    return ret;
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
    SparseSetT dimEOm, dimEOp, dimOm;

    // Intermediate values
    std::unique_ptr<Compute_wAOm> c1;

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

    c1.reset(new Compute_wAOm(dimIp, paramsA, gcmA, rmO));
    SparseSetT &dimAOp(c1->dimAOp);
    SparseSetT &dimAOm(c1->dimAOm);
    EigenColVectorT &wAOm_e(c1->wAOm_e);

    HntrSpec const &hntrA(cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr);
    HntrSpec const &hntrO(cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr);

printf("hntrO: %dx%d\n", hntrO.im, hntrO.jm);
printf("hntrA: %dx%d\n", hntrA.im, hntrA.jm);

    dimXAm.set_sparse_extent(X == 'A' ? gcmA->nA() : gcmA->nE());
    dimIp.set_sparse_extent(rmO->ice_regridder->nI());

    // Ensure that operations don't change universes.  Sanity check.
    std::unique_ptr<ConstUniverse> const_universe;

    if (X == 'E') {
        RegridParams paramsO(paramsA);
        paramsO.scale = false;
        paramsO.correctA = false;

        // Get the universe for EOp / EOm
        SparseSetT dimEOp;
        const_universe.reset(new ConstUniverse({"dimIp"}, {&dimIp}));
        linear::Weighted_Eigen &EOpvIp(tmp.take(rmO->matrix_d("EvI", {&dimEOp, &dimIp}, paramsO)));

        XOpvIp = &EOpvIp;

        const_universe.reset();        // Check that dimIp hasn't changed

        // ------------------------------------

        // EOpvAOp: Repeat A values in E
        const_universe.reset(new ConstUniverse({"dimEOp", "dimAOp"}, {&dimEOp, &dimAOp}));
        std::unique_ptr<linear::Weighted_Eigen> EOpvAOp(
            rmO->matrix_d("EvA", {&dimEOp, &dimAOp}, paramsO));

        const_universe.reset();        // Check that dims didn't change

        // Convert to EOmvAOm by removing cells not in dimAOm
        // EOmvAOm: Repeat A values in E
        blitz::Array<double,1> EOmvAOms(dimAOm.dense_extent());
        for (auto ii=EOpvAOp->M->begin(); ii != EOpvAOp->M->end(); ++ii) {
            auto const iAOp_d = ii->index(1);
            auto const iAOp_s = dimAOp.to_sparse(iAOp_d);
            int iAOm_d;
            if (dimAOm.to_dense_ignore_missing(iAOp_s, iAOm_d)) {
                auto const iEOp_d = ii->index(1);
                auto const iEOp_s = dimEOp.to_sparse(iEOp_d);
                auto const iEOm_d = dimEOp.add_dense(iEOp_s);  // add to dimEOp
                EOmvAOm_tl.add({iEOm_d, iAOm_d}, ii->value()); // add to EOmvAOm
                EOmvAOms(iAOm_d) = 1. / EOpvAOp->Mw(iAOp_d);   // add to EOmvAOmw
            }
        }



for (int i=0; i < 15; ++i) printf("EOmvAOms(%d) = %g\n", i, EOmvAOms(i));
for (int i=0; i < 15; ++i) printf("wAOm_e(%d) = %g\n", i, wAOm_e(i));

        // wEOm_e
        tmp.take<EigenColVectorT>(wXOm_e,
            EigenColVectorT(*EOmvAOm->M * map_eigen_diagonal(EOmvAOms) * wAOm_e));

        // ---------------- Compute the main matrix
        // ---------------- XAmvIp = XAmvXOm * XOpvIp
        // Rename variables
        SparseSetT &dimEAm(dimXAm);

        blitz::Array<double,1> wXOm(to_blitz(*wXOm_e));
printf("wXOm.extent(0) = %d\n", wXOm.extent(0));
for (int i=0; i < 15; ++i) printf("    wXOm[%d] = %g\n", i, wXOm(i));

        reset_ptr(XAmvXOm, MakeDenseEigenT(    // TODO: Call this XAvXO, since it's the same 'm' or 'p'
            std::bind(&raw_EOvEA, _1,
                hntrO, hntrA,
                eq_rad, &dimAOm, gcmA, wXOm),
            {SparsifyTransform::TO_DENSE_IGNORE_MISSING, SparsifyTransform::ADD_DENSE},
            {&dimEOm, &dimEAm}, 'T').to_eigen());
printf("XAmvXOm 1: %ld   dimEOm=%ld   dimEAm=%ld\n", (long)XAmvXOm->nonZeros(), dimEOm.dense_extent(), dimEAm.dense_extent());
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
    std::array<SparseSetT *,2> dimAm_array, dimAp_array;
    for (int i=0; i<2; ++i) dimAm_array[i] = nullptr;
    for (int i=0; i<2; ++i) dimAp_array[i] = nullptr;    // --> Don't change
    dimAm_array[index_Apdim] = &dimAm;      // Change the given index
    dimAp_array[index_Apdim] = &dimAp;
    TupleListT<2> Am_tl;
    spcopy(
        accum::to_dense_ignore_missing(dimAm_array,
        accum::to_sparse(dimAp_array,
        accum::ref(Am_tl))),
        Ap);
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
printf("BEGIN compute_XAmvIp\n");
    std::string dimXA(strprintf("dim%cA", X));
    std::unique_ptr<linear::Weighted_Eigen> ret(new linear::Weighted_Eigen(dims, false));    // not conservative

    // Do the legwork, and fish out things from the results that we need
    ComputeXAmvIp_Helper hh(dims, paramsA, gcmA, X, eq_rad, rmO);
    ret->tmp = std::move(hh.tmp);
    auto &XOpvIp(hh.XOpvIp);
    auto &wXAm_e(*hh.wXAm_e);
    auto &XAmvXOm(hh.XAmvXOm);

printf("XOpvIp: %ld\n", (long)(XOpvIp->M->nonZeros()));
printf("XAmvXOm: %ld\n", (long)(XAmvXOm->nonZeros()));

    // ----------- Put it all together (XAmvIp)
printf("BB1\n");
    blitz::Array<double,1> sXOpvIp(1. / XOpvIp->wM);
printf("BB2 |sXOpvIp|=%d %p\n", sXOpvIp.extent(0), sXOpvIp.data());
    ret->wM.reference(to_blitz(wXAm_e));    // wAAm_e or wEAm_e
printf("BB3\n");
    if (paramsA.scale) {
printf("AA1\n");
        auto sXAm(sum(*XAmvXOm, 0, '-'));
printf("AA2\n");
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(sXAm) *  *XAmvXOm *
            crop_mvp(hh.c1->dimAOm, hh.c1->dimAOp, 0,
                map_eigen_diagonal(sXOpvIp) * *XOpvIp->M)));
printf("AA3\n");
    } else {
printf("CC1 %p %p %p\n", &*ret, &*XAmvXOm, &*XOpvIp->M);
        ret->M.reset(new EigenSparseMatrixT(
            *XAmvXOm *
            crop_mvp(hh.c1->dimAOm, hh.c1->dimAOp, 0,
                map_eigen_diagonal(sXOpvIp) * *XOpvIp->M)));
printf("CC4\n");
    }
printf("BB4\n");
    ret->Mw.reference(XOpvIp->Mw);

printf("END compute_XAmvIp\n");
    return ret;
}
// -----------------------------------------------------------



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


    std::string const &matname(X=='A' ? "IvA" : "IvE");
    SparseSetT *dimIp(dims[0]);
    RegridParams paramsO(paramsA);
        paramsO.scale = false;
        paramsO.correctA = false;
    auto IpvXOp(
        rmO->matrix_d(matname,
            {dimIp, X == 'A' ? &hh.c1->dimAOp : &hh.dimEOp}, paramsO));

    // ----------- Put it all together (XAmvIp)
    blitz::Array<double,1> sXOpvIp(1. / XOpvIp->wM);

    auto XOmvXAm(XAmvXOm->transpose());
    auto &sXOmvXAm(hh.XAmvXOms);
    blitz::Array<double,1> sIpvXOp(1. / IpvXOp->wM);
    ret->wM.reference(XOpvIp->Mw);
    if (paramsA.scale) {
        ret->M.reset(new EigenSparseMatrixT(
            crop_mvp(hh.c1->dimAOm, hh.c1->dimAOp, 1,
                map_eigen_diagonal(sIpvXOp) * *IpvXOp->M) *
            map_eigen_diagonal(sXOmvXAm) * XOmvXAm
        ));
    } else {
        ret->M.reset(new EigenSparseMatrixT(
            crop_mvp(hh.c1->dimAOm, hh.c1->dimAOp, 1, *IpvXOp->M)
            * map_eigen_diagonal(sXOmvXAm) * XOmvXAm
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
    std::shared_ptr<icebin::GCMRegridder> const &_gcmO)
    :  gcmO(_gcmO)
{
    // Initialize superclass member
    agridA = make_agridA(gcmO->agridA);

    _ice_regridders = &gcmO->ice_regridders();

    correctA = gcmO->correctA;    // Not used
    hcdefs = gcmO->hcdefs;
    indexingHC = Indexing(
        {"O", "HC"},
        {0L, 0L},
        {agridA.dim.sparse_extent(), gcmO->indexingHC[1].extent},
        {gcmO->indexingHC.indices()[0], gcmO->indexingHC.indices()[1]});
    indexingE = derive_indexingE(agridA.indexing, indexingHC);

    auto const nO = gcmO->nA();
    foceanAOp.reference(blitz::Array<double,1>(nO));
    foceanAOm.reference(blitz::Array<double,1>(nO));
}



std::unique_ptr<RegridMatrices_Dynamic> GCMRegridder_ModelE::regrid_matrices(
    int sheet_index,
    blitz::Array<double,1> const &elevmaskI,
    RegridParams const &params) const
{
//    AbbrGrid const &gridO(gcmO->agridA);
    GridSpec_LonLat const &specO(cast_GridSpec_LonLat(*gcmO->agridA.spec));

    // Delicate construction of rmO
    std::unique_ptr<RegridMatrices_Dynamic> _rmO(
        gcmO->regrid_matrices(sheet_index, elevmaskI));
    std::unique_ptr<RegridMatrices_Dynamic> rm(
        new RegridMatrices_Dynamic(_rmO->ice_regridder, params));
    auto &rmO(rm->tmp.take(std::move(*_rmO)));

    rm->add_regrid("EAmvIp", std::bind(&compute_XAmvIp, _1, _2,
        this, 'E', specO.eq_rad, &rmO));
    rm->add_regrid("AAmvIp", std::bind(&compute_XAmvIp, _1, _2,
        this, 'A', specO.eq_rad, &rmO));


    rm->add_regrid("AAmvEAm", std::bind(*compute_AAmvEAm, _1, _2,
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


    rm->add_regrid("AvE", std::bind(*compute_AAmvEAm, _1, _2,
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
