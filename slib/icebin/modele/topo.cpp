#include <ibmisc/linear/eigen.hpp>
#include <icebin/ElevMask.hpp>
#include <icebin/modele/topo.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/grids.hpp>

using namespace blitz;
using namespace ibmisc;
using namespace spsparse;
using namespace std::placeholders;

static double const NaN = std::numeric_limits<double>::quiet_NaN();

namespace icebin {
namespace modele {



// ======================================================================
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


EigenColVectorT compute_wAOm(
    blitz::Array<double,1> const &foceanAOp,    // gcmA->foceanAOp
    blitz::Array<double,1> const &foceanAOm,    // gcmA->foceanAOm
    blitz::Array<double,1> const &wAOp,
    SparseSetT &dimAOp,   // Constant
    SparseSetT &dimAOm)   // write here
{
    ConstUniverseT const_dimAOp({"dimAOp"}, {&dimAOp});

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
// -------------------------------------------------------
// -------------------------------------------------------------
/** Helper class for raw_EOvEA().
@see raw_EOvEA */
class RawEOvEA {
public:
    MakeDenseEigenT::AccumT ret;    // Final place for EOvEA
    blitz::Array<double,1> const &wEO_d;
    // Things obtained from gcmA
    unsigned int const nhc;    // gcmA->nhc()
    ibmisc::Indexing const indexingHCO;    // gcmA->gcmO->indexingHC
    ibmisc::Indexing const indexingHCA;    // gcmA->indexingHC

int n=0;

    RawEOvEA(
        MakeDenseEigenT::AccumT &&_ret,
        blitz::Array<double,1> const &_wEO_d,
        unsigned int const _nhc,    // gcmA->nhc()
        ibmisc::Indexing const _indexingHCO,    // gcmA->gcmO->indexingHC
        ibmisc::Indexing const _indexingHCA)    // gcmA->indexingHC
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
void raw_EOvEA(
MakeDenseEigenT::AccumT &&ret,        // {dimEA, dimEO}; dimEO should not change here.
HntrSpec const &hntrO,
HntrSpec const &hntrA,
double const eq_rad,
SparseSetT const *dimAO,            // Used to clip in Hntr::matrix()
blitz::Array<double,1> &wEO_d,            // == EOvI.wM.  Dense indexing.
// Things obtained from gcmA
unsigned int const nhc,    // gcmA->nhc()
ibmisc::Indexing const indexingHCO,    // gcmA->gcmO->indexingHC
ibmisc::Indexing const indexingHCA)    // gcmA->indexingHC
{
    // Call Hntr to generate AOvAA; and use that (above) to produce EOvEA
    Hntr hntr_AOvAA(17.17, hntrO, hntrA, 0);    // dimB=A,  dimA=O

    hntr_AOvAA.overlap<RawEOvEA, DimClip>(
        RawEOvEA(std::move(ret), wEO_d, nhc, indexingHCO, indexingHCA),
        eq_rad, DimClip(dimAO));
}


// ------------------------------------------------------------------------
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
    ConstUniverseT const_dimAOm({"dimAOm"}, {&dimAOm});

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
// ------------------------------------------------------------------------

std::unique_ptr<linear::Weighted_Eigen> _compute_AAmvEAm(
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
    blitz::Array<double,1> const &wAOp)
{
    unsigned long const nA = indexingHCA[0].extent; // hntrA.size()
    unsigned int const nhc = indexingHCA[1].extent;
    unsigned long const nE = indexingHCA.extent();

    ConstUniverseT const_dimAOp({"dimEOp", "dimAOp"}, {&dimEOp, &dimAOp});

    SparseSetT &dimAAm(*dims[0]);
    SparseSetT &dimEAm(*dims[1]);

    // Must set sparse_extent
    dimAAm.set_sparse_extent(nA);
    dimEAm.set_sparse_extent(nE);

    std::unique_ptr<linear::Weighted_Eigen> ret(new linear::Weighted_Eigen(dims, false));    // not conservative

    // Compute wAOm (from wAOp)
    SparseSetT dimAOm;
    EigenColVectorT wAOm_e(compute_wAOm(foceanAOp, foceanAOm, wAOp, dimAOp, dimAOm));

    // ------------- Compute wAAm and AAmvAOm
    // Obtain hntrA and hntrO for process below.
    if (!hntrA.is_set()) (*icebin_error)(-1, "hntrA must be set");
    if (!hntrO.is_set()) (*icebin_error)(-1, "hntrO must be set");

    Hntr hntr_AOmvAAm(17.17, hntrO, hntrA);
    EigenSparseMatrixT AAmvAOm(MakeDenseEigenT(
        std::bind(&Hntr::overlap<MakeDenseEigenT::AccumT,DimClip>,
            &hntr_AOmvAAm, _1, eq_rad, DimClip(&dimAOm)),
        {SparsifyTransform::TO_DENSE_IGNORE_MISSING, SparsifyTransform::ADD_DENSE},
        {&dimAOm, &dimAAm}, 'T').to_eigen());

    blitz::Array<double,1> AAmvAOms(sum(AAmvAOm, 1, '-'));    // Note unusual way we weight/scale here
    auto &wAAm_e(ret->tmp.make<EigenColVectorT>(
        AAmvAOm * map_eigen_diagonal(AAmvAOms) * wAOm_e));

    // ------------ Compute AOmvEOm (from EOpvAOp)
    SparseSetT dimEOm;
    EigenSparseMatrixT EOmvAOm(compute_EOmvAOm_unscaled(dimEOm, dimAOm, EOpvAOp, dimEOp, dimAOp));
    auto EOmvAOms(sum(EOmvAOm, 1, '-'));
    auto AOmvEOm(EOmvAOm.transpose());
    auto &sAOmvEOm(EOmvAOms);

    // ------------- Compute wEOm
    EigenColVectorT wEOm_e(EOmvAOm * map_eigen_diagonal(EOmvAOms) * wAOm_e);
    auto wEOm(to_blitz(wEOm_e));

    // ------------ Compute EOmvEAm
    EigenSparseMatrixT EOmvEAm(MakeDenseEigenT(    // TODO: Call this EAvEO, since it's the same 'm' or 'p'
        std::bind(&raw_EOvEA, _1,
            hntrO, hntrA,
            eq_rad, &dimAOm, wEOm,
            nhc, indexingHCO, indexingHCA),
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
    if (scale) {
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(sAAmvAOm) * AAmvAOm *    // Works on whole cells, not just ice sheet; so we need to normalize by sAAmvAOm, not sAAm
            map_eigen_diagonal(sAOmvEOm) * AOmvEOm *
            map_eigen_diagonal(sEOmvEAm) * EOmvEAm));
    } else {
        blitz::Array<double,1> scale(ret->wM * sAAmvAOm);
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(scale) *
            AAmvAOm *    // Works on whole cells, not just ice sheet; so we need to normalize by sAAmvAOm, not sAAm
            map_eigen_diagonal(sAOmvEOm) * AOmvEOm *
            map_eigen_diagonal(sEOmvEAm) * EOmvEAm));
    }
    ret->Mw.reference(to_blitz(wEAm_e));

    return ret;
}

void make_topoA(
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
std::vector<uint16_t> const &underice_hc,    // gcmA->underice
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
blitz::Array<uint16_t,3> &underice3)
{
    ConstUniverseT const_dimAOp({"dimEOp", "dimAOp"}, {&dimEOp, &dimAOp});


    Hntr hntr_AvO(17.17, hspecA, hspecO);

    blitz::Array<double, 2> WTO(const_array(blitz::shape(hspecO.jm,hspecO.im), 1.0));
    hntr_AvO.regrid(WTO, foceanOm2, foceanA2);
    hntr_AvO.regrid(WTO, flakeOm2, flakeA2);
    hntr_AvO.regrid(WTO, fgrndOm2, fgrndA2);
    hntr_AvO.regrid(WTO, fgiceOm2, fgiceA2);
    hntr_AvO.regrid(WTO, zatmoOm2, zatmoA2);
    hntr_AvO.regrid(WTO, zlakeOm2, zlakeA2);
    hntr_AvO.regrid(fgiceOm2, zicetopOm2, zicetopA2);

    auto foceanOp(reshape1(foceanOp2));
    auto foceanOm(reshape1(foceanOm2));
    auto flakeOm(reshape1(flakeOm2));
    auto fgiceOm(reshape1(fgiceOm2));
    auto zatmoOm(reshape1(zatmoOm2));
    auto zicetopO(reshape1(zicetopOm2));

    auto foceanA(reshape1(foceanA2));
    auto flakeA(reshape1(flakeA2));
    auto fgrndA(reshape1(fgrndA2));
    auto fgiceA(reshape1(fgiceA2));
    auto zatmoA(reshape1(zatmoA2));

    // ================= Create fhc, elevE and underice
    int const nhc_icebin = hcdefs.size();
    int const nhc_gcm = get_nhc_gcm(nhc_icebin);
    auto const nA = indexingHCA[0].extent;
    blitz::TinyVector<int,2> shapeE2(nhc_gcm, nA);

    // Initialize
    fhc3 = 0;
    elevE3 = NaN;
    underice3 = 0;    // Elevation grid cell unused

    // ------------ Segment 0: Elevation Classes
    printf("Segment 0: Elevation Classes\n");
    blitz::Array<double,2> fhcE2(reshape<double,3,2>(fhc3, shapeE2));
    auto elevE2(reshape<double,3,2>(elevE3, shapeE2));
    auto undericeE2(reshape<uint16_t,3,2>(underice3, shapeE2));

    std::array<int,2> iTuple;
        int &iA2(iTuple[0]);
        int &ihc(iTuple[1]);

    // Compute AOmvEOm --> fhc
    auto wAOp(sum(EOpvAOp, 1, '+'));
    SparseSetT dimAAm,dimEAm;
    std::unique_ptr<linear::Weighted_Eigen> AAmvEAm(_compute_AAmvEAm(
        {&dimAAm, &dimEAm}, true, eq_rad,    // scale=true
        hspecO, hspecA, indexingHCO, indexingHCA,
        foceanOp, foceanOm,
        EOpvAOp, dimEOp, dimAOp, wAOp));

    for (auto ii=begin(*AAmvEAm->M); ii != end(*AAmvEAm->M); ++ii) {
        int const iA_d = ii->index(0);
        int const iA = dimAAm.to_sparse(iA_d);
        int const iE_d = ii->index(1);
        int const iE = dimEAm.to_sparse(iE_d);

        indexingHCA.index_to_tuple(&iTuple[0], iE);

        // iE must be contained within cell iA (local propety of matrix)
        if (iA2 != iA) (*icebin_error)(-1,
            "Matrix is non-local: iA=%ld, iE=%ld, iA2=%ld",
            (long)iA, (long)iE, (long)iA2);

        if (ihc < 0 || ihc >= nhc_icebin) (*icebin_error)(-1,
            "ihc out of range [0,%d): %d", nhc_icebin, ihc);

        if (iA < 0 || iA >= shapeE2(1)) (*icebin_error)(-1,
            "iA out of range [0,%d): %d", shapeE2(1), iA);

//if (values(i) < 0 || values(i) >= 1.) printf("AvE(%d, ihc=%d) = %g\n", iA, ihc, values(i));

        int const _ihc = ihc;    // Get around bug in Blitz++
        fhcE2(ihc,iA) += ii->value();
        undericeE2(ihc,iA) = underice_hc[ihc];    // No IceBin coupling here
        elevE2(ihc, iA) = hcdefs[ihc];
    }

    // ------------ Segment 1: land part of sealand
    printf("Segment 1: seaLAND\n");
    int ec_base = nhc_icebin;    
    for (int j=0; j<hspecA.jm; ++j) {
    for (int i=0; i<hspecA.im; ++i) {
        if (fgiceA(j,i) > 0) {
            fhc3(ec_base, j,i) = (fgiceA2(j,i) == 0 ? 0 : 1e-30);
            elevE3(ec_base, j,i) = zicetopA2(j,i);
            underice3(ec_base, j,i) = UI_NOTHING;
        }
    }}
}


}}    // namespace icebin::modele
