#include <ibmisc/linear/eigen.hpp>
#include <ibmisc/stdio.hpp>
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

    // It would be more correct to write wAOp * AOmvAOp
    // But scaled_AOmvAOp is diagonal, so it doesn't matter here.
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
std::unique_ptr<linear::Weighted_Eigen> _compute_AAmvEAm_EIGEN(
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
    EigenSparseMatrixT const &EOpvAOp,    // unscaled
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

    std::unique_ptr<linear::Weighted_Eigen> AAmvEAm(
        new linear::Weighted_Eigen({&dimAAm,&dimEAm}, false));    // not conservative

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
    auto &wAAm_e(AAmvEAm->tmp.make<EigenColVectorT>(
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
    EigenSparseMatrixT EOmvEAm(MakeDenseEigenT(    // TODO: Call this EOvEA, since it's the same 'm' or 'p'
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
    auto &wEAm_e(AAmvEAm->tmp.make<EigenColVectorT>(
        EAmvEOm * map_eigen_diagonal(EAmvEOms) * wEOm_e));


    // ------------- Put it all together
    AAmvEAm->wM.reference(to_blitz(wAAm_e));
    blitz::Array<double,1> sAAmvAOm(sum(AAmvAOm, 0, '-'));
    if (scale) {
        AAmvEAm->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(sAAmvAOm) * AAmvAOm *    // Works on whole cells, not just ice sheet; so we need to normalize by sAAmvAOm, not sAAm
            map_eigen_diagonal(sAOmvEOm) * AOmvEOm *
            map_eigen_diagonal(sEOmvEAm) * EOmvEAm));
    } else {
        // Maintains identity: M_unscaled = wM * M_scaled
        blitz::Array<double,1> scale(AAmvEAm->wM * sAAmvAOm);
        AAmvEAm->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(scale) *
            AAmvAOm *    // Works on whole cells, not just ice sheet; so we need to normalize by sAAmvAOm, not sAAm
            map_eigen_diagonal(sAOmvEOm) * AOmvEOm *
            map_eigen_diagonal(sEOmvEAm) * EOmvEAm));
    }
    AAmvEAm->Mw.reference(to_blitz(wEAm_e));

    return AAmvEAm;

}

linear::Weighted_Tuple _compute_AAmvEAm(
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
    EigenSparseMatrixT const &EOpvAOp,    // unscaled
    SparseSetT &dimEOp,
    SparseSetT &dimAOp,
    blitz::Array<double,1> const &wAOp)
{
    SparseSetT dimAAm, dimEAm;
    std::unique_ptr<linear::Weighted_Eigen> ret(_compute_AAmvEAm_EIGEN(
        {&dimAAm, &dimEAm},
        scale, eq_rad, hntrO, hntrA,
        indexingHCO, indexingHCA, foceanAOp, foceanAOm,
        EOpvAOp, dimEOp, dimAOp, wAOp));
    return to_tuple(*ret);
}

/** Create, allocate and load data into a bundle representing a TOPOO file.
@param type Allows for variations on which variables are added to the bundle
    BundleOType::MERGEO: Variables appropriate for make_merged_topoo.cpp
    BundleOType::MAKEA: Variables appropriate for input of make_topoa.cpp
@param topoO_fname
    Name of TOPOO file to load (output of make_topoo.cpp or make_merged_topoo.cpp)
    If not set, then don't load or allocate anything.
*/
ibmisc::ArrayBundle<double,2> topoo_bundle(
BundleOType type,
std::string const &topoO_fname)
{
    ibmisc::ArrayBundle<double,2> topoo;

    // ------------- Non-rounded versions (Op)
    topoo.add("FOCEANF", {
        "description", "Fractional ocean ocver",
        "units", "1",
        "sources", "GISS 1Qx1",
    });
    if (type == BundleOType::MERGEO) {
        auto &fgiceOp(topoo.add("FGICEF", {
            "description", "Glacial Ice Surface Fraction (Ocean NOT rounded)",
            "units", "0:1",
            "sources", "GISS 1Qx1",
        }));
        auto &zatmoOp(topoo.add("ZATMOF", {
            "description", "Atmospheric Topography",
            "units", "m",
            "sources", "ETOPO2 1Qx1",
        }));
    }



    // ------------ Rounded Versions (Om)
    topoo.add("FOCEAN", {
        "description", "0 or 1, Bering Strait 1 cell wide",
        "units", "1",
        "source", "GISS 1Qx1",
    });

    auto &flakeOm(topoo.add("FLAKE", {
        "description", "Lake Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &fgrndOm(topoo.add("FGRND", {
        "description", "Ground Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &fgiceOm(topoo.add("FGICE", {
        "description", "Glacial Ice Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &zatmoOm(topoo.add("ZATMO", {
        "description", "Atmospheric Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &zlakeOm(topoo.add("ZLAKE", {
        "description", "Lake Surface Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &zicetopOm(topoo.add("ZICETOP", {
        "description", "Atmospheric Topography (Ice-Covered Regions Only)",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));

//    // ZLAND_MIN/ZLAND_MAX are not computed for global ice
//    if (type != BundleOType::MERGEO) {
        auto &zland_minOm(topoo.add("ZLAND_MIN", {
            "description", "Minimum land/ice topography in this cell",
            "units", "m",
            "sources", "ETOPO2 1Qx1",
        }));
        auto &zland_maxOm(topoo.add("ZLAND_MAX", {
            "description", "Maximum land/ice topography in this cell",
            "units", "m",
            "sources", "ETOPO2 1Qx1",
        }));
//    }

    // Read TOPOO input
    if (topoO_fname != "") {
        NcIO topoo_nc(topoO_fname, 'r');

        // Read from topoO file, and allocate resulting arrays.
        topoo.ncio_alloc(topoo_nc, {}, "", "double",
            get_dims(topoo_nc, {"jm", "im"}));
    }

    return topoo;

}
// ------------------------------------------------------------------------
TopoABundles::TopoABundles(
ibmisc::ArrayBundle<double,2> const &topoo,
HntrSpec const &hspecA,
int const nhc_gcm)
{

    static std::vector<std::string> const varsO {"FOCEAN", "FLAKE", "FGRND", "FGICE", "ZATMO", "ZLAKE", "ZICETOP", "ZLAND_MIN", "ZLAND_MAX"};
    static std::vector<std::string> const varsA {"focean", "flake", "fgrnd", "fgice", "zatmo", "hlake", "zicetop", "zland_min", "zland_max"};

    for (size_t i=0; i<varsO.size(); ++i) {
        std::string const &nameO(varsO[i]);
        std::string const &nameA(varsA[i]);

        // Get data record for this variable
        auto &d(topoo.at(nameO));

        // Add to topoa based on info from topoo
        this->a.add(ibmisc::ArrayBundle<double,2>::Data(
            nameA, blitz::Array<double,2>(),
            std::array<int,2>{hspecA.jm, hspecA.im},
            d.meta.sdims,
            std::vector<std::pair<std::string, std::string>>(d.meta.attr)));
    }
    this->a_i.add("mergemask", {
        "description", "1 where TOPO generator has merged from local ice sheets"
    });

    this->a.allocate(std::array<int,2>{hspecA.jm, hspecA.im}, {"jm","im"}, false);    // check=false
    this->a_i.allocate(std::array<int,2>{hspecA.jm, hspecA.im}, {"jm","im"}, false);    // check=false

    // --------------- Allocate 3D arrays to go in TOPOA file
    auto &fhc(this->a3.add("fhc", {
        "description", "fraction of ice-covered area for each elevation class",
        "units", "1"
    }));
    auto &elevE(this->a3.add("elevE", {
        "description", "Elevation of each elevation class",
        "units", "1"
    }));
    auto &underice(this->a3_i.add("underice", {
        "description", "Model below the show/firn (UI_UNUSED=0, UI_LOCALICE=1, UI_GLOBALICE=2, UI_VGHOST=3, UI_HGHOST=4, UI_SEALAND=4)"
    }));

    std::array<int,3> shape3 {nhc_gcm, hspecA.jm, hspecA.im};
    this->a3.allocate(shape3, {"nhc", "jm", "im"});
    this->a3_i.allocate(shape3, {"nhc", "jm", "im"});

}

void TopoABundles::ncio(NcIO &ncio, std::vector<netCDF::NcDim> const &ncdims)
{
    auto &nhc(ncdims[0]);
    auto &jm(ncdims[1]);
    auto &im(ncdims[2]);

    this->a.ncio(ncio, {}, "", "double", {jm, im});
    this->a_i.ncio(ncio, {}, "", "short", {jm,im});  // Must be short for NetCDF3
    this->a3.ncio(ncio, {}, "", "double", ncdims);
    this->a3_i.ncio(ncio, {}, "", "short", ncdims);  // Must be short for NetCDF3
}

// ------------------------------------------------------------------------
static void merge_poles(blitz::Array<double,2> &var)
{
    std::array<int,2> const jix {var.lbound(0), var.ubound(0)};
    for (int j : jix) {
        double sum = 0;
        for (int i=var.lbound(1); i <= var.ubound(1); ++i)
            sum += var(j,i);
        double const mean = sum / (double)(var.ubound(1)-var.lbound(1)+1);

        for (int i=var.lbound(1); i <= var.ubound(1); ++i)
            var(j,i) = mean;
    }
}
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// This is an spsparse accumulator, gets called on every
// matrix element by hntr_AvO.scaled_regrid_matrix() below.
struct _RegridMinMax {
    // See spsparse/spsparse.hpp AccumTraits
    static int const rank = 2;
    typedef int index_type;
    typedef double val_type;

    blitz::Array<double,1> zland_minO, zland_maxO;
    blitz::Array<double,1> zland_minA, zland_maxA;
    blitz::Array<int16_t,1> mergemaskO, mergemaskA;

    void add(std::array<int,2> const &index, double val)
    {
        int const iA = index[0];
        int const iO = index[1];

        // mergemaskA=1 if any of associated mergemaskO is 1
        if (mergemaskO(iO)) {
            mergemaskA(iA) = 1;

            // minA = min(all minO)
            zland_minA(iA) = std::min(zland_minA(iA), zland_minO(iO));
            zland_maxA(iA) = std::max(zland_maxA(iA), zland_maxO(iO));
        }
    }
};
// ------------------------------------------------------------------------
std::vector<std::string> make_topoA(
// AAmvEAM is either read from output of global_ec (for just global ice);
// or it's the output of compute_AAmvEAm_merged (for merged global+local ice)
blitz::Array<double,2> const &foceanOm2,     // Rounded FOCEAN
blitz::Array<double,2> const &flakeOm2,
blitz::Array<double,2> const &fgrndOm2,
blitz::Array<double,2> const &fgiceOm2,
blitz::Array<double,2> const &zatmoOm2,
blitz::Array<double,2> const &zlakeOm2,
blitz::Array<double,2> const &zicetopOm2,
blitz::Array<double,2> const &zland_minOm2,
blitz::Array<double,2> const &zland_maxOm2,
blitz::Array<int16_t,2> const &mergemaskOm2,
//
// Things obtained from gcmA
HntrSpec const &hspecO,        // cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr
HntrSpec const &hspecA,        // cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr
ibmisc::Indexing const indexingHCA,    // gcmA->indexingHC
std::vector<double> const &hcdefs,        // gcmA->hcdefs()
//int const nhc_icesheet,                         // gcmA->gcmO->hcdefs()Just EC's related to dynamic ice sheets
std::vector<int16_t> const &underice_hc,    // gcmA->underice
//
linear::Weighted_Tuple const &AAmvEAm,
// ----- Outputs; Native NetCDF indexing (j,i) zero-based
blitz::Array<double,2> &foceanA2,    // Rounded FOCEAN
blitz::Array<double,2> &flakeA2,
blitz::Array<double,2> &fgrndA2,
blitz::Array<double,2> &fgiceA2,
blitz::Array<double,2> &zatmoA2,
blitz::Array<double,2> &zlakeA2,
blitz::Array<double,2> &zicetopA2,
blitz::Array<double,2> &zland_minA2,
blitz::Array<double,2> &zland_maxA2,
blitz::Array<int16_t,2> &mergemaskA2,
//
blitz::Array<double,3> &fhc3,
blitz::Array<double,3> &elevE3,
blitz::Array<int16_t,3> &underice3)
{

    Hntr hntr_AvO(17.17, hspecA, hspecO);

    blitz::Array<double, 2> WTO(const_array(blitz::shape(hspecO.jm,hspecO.im), 1.0));
    hntr_AvO.regrid(WTO, foceanOm2, foceanA2);
    hntr_AvO.regrid(WTO, flakeOm2, flakeA2);
    hntr_AvO.regrid(WTO, fgrndOm2, fgrndA2);
    hntr_AvO.regrid(WTO, fgiceOm2, fgiceA2);
    hntr_AvO.regrid(WTO, zatmoOm2, zatmoA2);
    hntr_AvO.regrid(WTO, zlakeOm2, zlakeA2);
    hntr_AvO.regrid(fgiceOm2, zicetopOm2, zicetopA2);

    // -------------------------
    // Regrid mergemask (mask, not a double)
    // Also set zland_min and zland_max

    // These vars aren't set everywhere in regridding, so they must be
    // initialized.
    mergemaskA2 = 0;
    zland_minA2 = std::numeric_limits<double>::max();
    zland_maxA2 = std::numeric_limits<double>::min();
    {_RegridMinMax rmm;
        rmm.zland_minO.reference(reshape1(zland_minOm2));
        rmm.zland_maxO.reference(reshape1(zland_maxOm2));
        rmm.zland_minA.reference(reshape1(zland_minA2));
        rmm.zland_maxA.reference(reshape1(zland_maxA2));
        rmm.mergemaskO.reference(reshape1(mergemaskOm2));
        rmm.mergemaskA.reference(reshape1(mergemaskA2));

        hntr_AvO.scaled_regrid_matrix(accum::ref(rmm));
    }
    for (int j=0; j<hspecA.jm; ++j) {
    for (int i=0; i<hspecA.im; ++i) {
        if (zland_minA2(j,i) == std::numeric_limits<double>::max()) zland_minA2(j,i) = NaN;
        if (zland_maxA2(j,i) == std::numeric_limits<double>::min()) zland_maxA2(j,i) = NaN;
    }}
    // -------------------------

    merge_poles(foceanA2);
    merge_poles(flakeA2);
    merge_poles(fgrndA2);
    merge_poles(fgiceA2);
    merge_poles(zatmoA2);
    merge_poles(zlakeA2);
    merge_poles(zicetopA2);
    merge_poles(zland_minA2);
    merge_poles(zland_maxA2);


#if 0    // not needed
    auto foceanOm(reshape1(foceanOm2));
    auto flakeOm(reshape1(flakeOm2));
    auto fgiceOm(reshape1(fgiceOm2));
    auto zatmoOm(reshape1(zatmoOm2));
    auto zicetopO(reshape1(zicetopOm2));
#endif

    auto foceanA(reshape1(foceanA2));
    auto flakeA(reshape1(flakeA2));
    auto fgrndA(reshape1(fgrndA2));
    auto fgiceA(reshape1(fgiceA2));
    auto zatmoA(reshape1(zatmoA2));

    // ================= Create fhc, elevE and underice
    int const nhc_icebin = hcdefs.size();   // local + global EC's
    int nhc_local;    // Number of ECs devoted to local ice sheets (not global ice)
    for (nhc_local=0; nhc_local<nhc_icebin; ++nhc_local) {
        // underice_hc contains UI_LOCALICE for local ice; and UI_GLOBALICE for global ice.
        // See merge_topo.cpp for how underice_hc is constructed.
        if (underice_hc[nhc_local] != UI_LOCALICE) break;
    }
    int const nhc_gcm = get_nhc_gcm(nhc_icebin);
    auto const nA = indexingHCA[0].extent;
    blitz::TinyVector<int,2> shapeE2(nhc_gcm, nA);

    // Initialize
    fhc3 = 0;
    elevE3 = NaN;
    underice3 = UI_UNUSED;    // Elevation grid cell unused

    // ------------ Segment 0: Elevation Classes
    printf("Segment 0: Elevation Classes\n");
    blitz::Array<double,2> fhcE2(reshape<double,3,2>(fhc3, shapeE2));
    auto elevE2(reshape<double,3,2>(elevE3, shapeE2));
    auto undericeE2(reshape<int16_t,3,2>(underice3, shapeE2));

    {std::array<int,2> iTuple;
        int &iA2(iTuple[0]);
        int &ihc(iTuple[1]);

        for (auto ii=begin(AAmvEAm.M); ii != end(AAmvEAm.M); ++ii) {
            int const iA = ii->index(0);
            int const iE = ii->index(1);

            indexingHCA.index_to_tuple(&iTuple[0], iE);

            // iE must be contained within cell iA (local propety of matrix)
            if (iA2 != iA) (*icebin_error)(-1,
                "Matrix is non-local: iA=%ld, iE=%ld, iA2=%ld",
                (long)iA, (long)iE, (long)iA2);

            if (ihc < 0 || ihc >= nhc_icebin) (*icebin_error)(-1,
                "ihc out of range [0,%d): %d", nhc_icebin, ihc);

            if (iA < 0 || iA >= shapeE2(1)) (*icebin_error)(-1,
                "iA out of range [0,%d): %d", shapeE2(1), iA);

            int const _ihc = ihc;    // Get around bug in Blitz++
            fhcE2(ihc,iA) += ii->value();
            undericeE2(ihc,iA) = underice_hc[ihc];    // No IceBin coupling here
        }
    }

    // Merge grid cells on south pole
    // (there is no ice on north pole, so fhc==0 alrady there)
    for (int ihc=0; ihc<nhc_icebin; ++ihc) {
        double fhc_sum = 0;
        for (int i=0; i<hspecA.im; ++i) fhc_sum += fhc3(ihc,0,i);
        double const fhc_mean = fhc_sum / (double)hspecA.im;
        for (int i=0; i<hspecA.im; ++i) fhc3(ihc,0,i) = fhc_mean;
    }

printf("nhc_icebin = %d %d\n", nhc_icebin, nhc_local);
    // Set EC levels and EC halos (ghost ECs)
    for (int j=0; j<hspecA.jm; ++j) {
    for (int i=0; i<hspecA.im; ++i) {
        const int minhc0 = std::numeric_limits<int>::max();
        const int maxhc0 = std::numeric_limits<int>::min();
        int minhc = minhc0;
        int maxhc = maxhc0;
        int zland_minhc=std::numeric_limits<int>::max();
        int zland_maxhc=std::numeric_limits<int>::min();
        for (int ihc=0; ihc<nhc_local; ++ihc) {
            // Set elevE everywhere.  This is required by ModelE
            // See downscale_temperature_li(), which doesn't
            // look at fhc.
            elevE3(ihc,j,i) = hcdefs[ihc];

            // zland_minhc will end up being the largest EC level < zland_min
            if (elevE3(ihc,j,i) < zland_minA2(j,i)) zland_minhc = ihc;
            // zland_maxhc will end up being the largest EC level <= zland_max
            if (elevE3(ihc,j,i) <= zland_maxA2(j,i)) zland_maxhc = ihc;

            // Determine min and max EC for this Atmosphere gridcell
            if (fhc3(ihc,j,i) != 0) {
                minhc = std::min(minhc,ihc);
                maxhc = std::max(maxhc,ihc);
            }
        }

        // Only put ghost points on grid cells containing land.
        // Ghost points over pure ocean will cause:
        //      >>  PBL: Ts out of range  <<
        // in ModelE:PBL_DRV.f
        if (std::abs(foceanA2(j,i) - 1.0) > 1.e-14) {
            int minghost = std::max(0, zland_minhc-1);
            int maxghost = std::min(zland_maxhc+2, nhc_local-1);

if (!std::isnan(zland_minA2(j,i))) printf("minmax %d,%d (%f %f) (%d %d) (%d %d)\n", j,i, zland_minA2(j,i), zland_maxA2(j,i), zland_minhc, zland_maxhc, minghost, maxghost);

             for (int ihc=minghost; ihc<=maxghost; ++ihc) {
                 if (fhc3(ihc,j,i)==0) {
                     fhc3(ihc,j,i) = 1.e-30;
                     underice3(ihc,j,i) = UI_VGHOST;
                 }
             }
        }

        // Set elevE for global ice
        for (int ihc=nhc_local; ihc<nhc_icebin; ++ihc) {
            elevE3(ihc,j,i) = hcdefs[ihc];
        }


// #if 0
//         // If this cell has no ice, set minhc and maxhc based on
//         // ei_land
//         if (minhc == minhc0 && !std::isnan(zland_minA2(j,i))) {
//                 minghost = std::max(0, zland_minhc-1);
//                 maxghost = std::min(zland_maxhc+2, nhc_local-1);
// 
// 
//         }
// 
// 
//         // Make ghost points on cells that MIGHT have ice someday,
//         // because they are intersected by a local ice grid.
// #endif
// //        if (mergemaskA2(j,i)) {
// 
//             int minghost, maxghost;
// //            if (maxhc < 0) {
//                 // If the gridcell has no ice in it, set up ghost points at
//                 // elevations that MIGHT see ice.  (Assume Z interpolation in
//                 // IceBin).
//                 minghost = std::max(0, zland_minhc-1);
//                 maxghost = std::min(zland_maxhc+2, nhc_local-1);
// //            } else {
// //                minghost = std::max(0,minhc-2);
// //                maxghost = std::min(maxhc+2,nhc_local-1);
// //            }
// 
//             for (int ihc=minghost; ihc<=maxghost; ++ihc) {
//                 if (fhc3(ihc,j,i)==0) {
//                     fhc3(ihc,j,i) = 1.e-30;
//                     underice3(ihc,j,i) = UI_VGHOST;
//                 }
//             }
// //        }
     }}


    // ------------ Segment 1: land part of sealand
    printf("Segment 1: seaLAND\n");
    int ec_base = nhc_icebin;    
    for (int j=0; j<hspecA.jm; ++j) {
    for (int i=0; i<hspecA.im; ++i) {
        if (fgiceA2(j,i) > 0) {
            fhc3(ec_base, j,i) = (fgiceA2(j,i) == 0 ? 0 : 1e-30);
            underice3(ec_base, j,i) = UI_SEALAND;
        }

        // elevE must be non-NaN everywhere, because ModelE
        // downscales temperature everywhere, for all ECs
        // See SURFACE_LANDICE.F90:downscale_temperature_li()
        elevE3(ec_base, j,i) = zatmoA2(j,i);
    }}

    // ==================================================
    // Sanity Check the Result
    std::vector<std::string> errors;
    sanity_check_land_fractions(foceanA2, flakeA2, fgrndA2, fgiceA2, errors);
    sanity_check_fhc(fhc3, errors);

    return errors;
}

void sanity_check_fhc(
blitz::Array<double,3> const &fhc3,
std::vector<std::string> &errors)
{
    for (int j=0; j<fhc3.extent(1); ++j) {
    for (int i=0; i<fhc3.extent(2); ++i) {
        // ---------- FHC must add up to 1
        double all_fhc = 0;
        for (int ihc=0; ihc<fhc3.extent(0); ++ihc) all_fhc += fhc3(ihc,j,i);
        all_fhc += 1.0;    // Disregard 1e-30 values
        if (all_fhc != 1.0 && std::abs(all_fhc-2.0) > 1.e-13) errors.push_back(
            ibmisc::strprintf("(%d, %d): sum(FHC) = %g", i+1,j+1, all_fhc-1.0));
    }}

}

void sanity_check_land_fractions(
blitz::Array<double,2> const &foceanA2,    // Rounded FOCEAN
blitz::Array<double,2> const &flakeA2,
blitz::Array<double,2> const &fgrndA2,
blitz::Array<double,2> const &fgiceA2,
std::vector<std::string> &errors)
{
    for (int j=0; j<foceanA2.extent(0); ++j) {
    for (int i=0; i<foceanA2.extent(1); ++i) {
        // ----------- Land fractions must add up to 1
        double const all_frac = foceanA2(j,i) + fgrndA2(j,i) + flakeA2(j,i) + fgiceA2(j,i);
        if (std::abs(all_frac-1.0) > 1.e-13) errors.push_back(
            ibmisc::strprintf("(%d, %d): FOCEAN(%g) + FGRND(%g) + FLAKE(%g) + FGICE(%g)  = %g",
            i+1,j+1, foceanA2(j,i), fgrndA2(j,i), flakeA2(j,i), fgiceA2(j,i), all_frac));
    }}
}


}}    // namespace icebin::modele
