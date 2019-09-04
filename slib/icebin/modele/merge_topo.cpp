#include <icebin/modele/merge_topo.hpp>
#include <icebin/modele/topo.hpp>
#include <icebin/modele/grids.hpp>
#include <icebin/eigen_types.hpp>
#include <ibmisc/linear/compressed.hpp>

using namespace blitz;
using namespace ibmisc;
using namespace spsparse;

static double const NaN = std::numeric_limits<double>::quiet_NaN();

namespace icebin {
namespace modele {

/** Return value for get_sheet_elevO() */
class GetSheetElevO {
public:
    TmpAlloc _tmp;
    SparseSetT dimO;    // Dense indexing used for wO and elevO
    blitz::Array<double,1> wO;        // Ice (or land) covered area of each gridcell [m^2]
    blitz::Array<double,1> elev_areaO;     // Elevation [m] * Ice-covered Cell Area [m^2] = [m^3]
};

/** Returns elevation on ocean grid (elevO) for a single local ice sheet.
@param gcmO GCMRegridder
@param paramsA Regrid parameters.
@param sheet_index Index of the ice sheet within gcmO for which we seek an answer.
@param elevmaskI Combined elevation/mask for the selected ice sheet.
    =elev or NaN; for either all land, or ice-covered land, depending
    on desired result. */
static GetSheetElevO get_sheet_elevO(
GCMRegridder_Standard *gcmO,
RegridParams const &paramsO,
int sheet_index,
blitz::Array<double,1> const &elevmaskI)
{
    GetSheetElevO ret;

    size_t const nI = gcmO->nI(sheet_index);
    size_t const nO = gcmO->nA();

    // Obtain the OvI matrix
    SparseSetT dimI;
    std::unique_ptr<RegridMatrices_Dynamic> rmO(
        gcmO->regrid_matrices(sheet_index, elevmaskI, paramsO));
   std::unique_ptr<ibmisc::linear::Weighted_Eigen> OvI(
        rmO->matrix_d("AvI", {&ret.dimO, &dimI}, paramsO));


    // Construct elevI in dense indexing space
    blitz::Array<double,1> elevI(dimI.dense_extent());
    for (size_t iI_d=0; iI_d < dimI.dense_extent(); ++iI_d) {
        auto iI_s = dimI.to_sparse(iI_d);
        elevI(iI_d) = elevmaskI(iI_s); //emI.elev(iI_s);
    }

    ret.elev_areaO.reference(OvI->apply(elevI, NaN, false, ret._tmp));    // dense indexing, force_conservation=false
    ret.wO.reference(OvI->wM);

    return ret;
}

void sanity_nonan(
std::string const &label,
blitz::Array<double,2> const &varA2,    // Rounded FOCEAN
std::vector<std::string> &errors)
{
    for (int j=0; j<varA2.extent(0); ++j) {
    for (int i=0; i<varA2.extent(1); ++i) {
        if (std::isnan(varA2(j,i))) errors.push_back(
            ibmisc::strprintf("(%d, %d): %s is NaN", i+1,j+1, label.c_str()));
    }}
}


void merge_topoO(
// ------ TOPOO arrays, originally with just global ice
// Ice sheets will be merged into them.
//
// Ice model viewpoint (fractional ocean cells)
blitz::Array<double,2> &foceanOp2,    // Fractional FOCEAN
blitz::Array<double,2> &fgiceOp2,
blitz::Array<double,2> &zatmoOp2,
// ModelE viewpoint (rounded ocean cells)
blitz::Array<double,2> &foceanOm2,     // Rounded FOCEAN
blitz::Array<double,2> &flakeOm2,
blitz::Array<double,2> &fgrndOm2,
blitz::Array<double,2> &fgiceOm2,
blitz::Array<double,2> &zatmoOm2,
// Not affected by Om; as long as top of ice is maintained even for ocean-rounded cells.
blitz::Array<double,2> &zicetopO2,
// ------ Local ice sheets to merge in...
GCMRegridder_Standard *gcmO,    // Multiple IceRegridders
RegridParams const &paramsA,
std::vector<blitz::Array<double,1>> const &emI_lands,
std::vector<blitz::Array<double,1>> const &emI_ices,
double const eq_rad,    // Radius of the earth
std::vector<std::string> &errors)
{

    auto foceanOp(reshape1(foceanOp2));
    auto fgiceOp(reshape1(fgiceOp2));
    auto zatmoOp(reshape1(zatmoOp2));
    auto foceanOm(reshape1(foceanOm2));
    auto flakeOm(reshape1(flakeOm2));
    auto fgrndOm(reshape1(fgrndOm2));
    auto fgiceOm(reshape1(fgiceOm2));
    auto zatmoOm(reshape1(zatmoOm2));
    auto zicetopO(reshape1(zicetopO2));

    // Sanity check inputs
    sanity_nonan("foceanOp2-0", foceanOp2, errors);
    sanity_nonan("fgiceOp2-0", fgiceOp2, errors);
    sanity_nonan("zatmoOp2-0", zatmoOp2, errors);
    sanity_nonan("foceanOm2-0", foceanOm2, errors);
    sanity_nonan("flakeOm2-0", flakeOm2, errors);
    sanity_nonan("fgrndOm2-0", fgrndOm2, errors);
    sanity_nonan("fgiceOm2-0", fgiceOm2, errors);
    sanity_nonan("zatmoOm2-0", zatmoOm2, errors);
    sanity_nonan("zicetopO2-0", zicetopO2, errors);

    auto nO = gcmO->nA();


    // Accumulate result of ice sheets
    blitz::Array<double,1> da_giceO(nO);     // increase in area of landice [m^2]
    da_giceO = 0.;
    blitz::Array<double,1> da_zicetopO(nO);     // increase in area of landice * elevation [m^3]
    da_zicetopO = 0.;
    blitz::Array<double,1> da_contO(nO);    // increase in area of continent [m^2]
    da_contO = 0.;
    blitz::Array<double,1> da_zatmoO(nO);    // increase in area of continent * elevation [m^3]
    da_zatmoO = 0.;

    RegridParams paramsO_rawA(paramsA);
        paramsO_rawA.scale = true;
        paramsO_rawA.correctA = false;
    RegridParams paramsO_correctA(paramsA);
        paramsO_correctA.scale = false;
        paramsO_correctA.correctA = true;

    // -------------- Compute additional area of land and ice due to local ice sheets
    for (size_t sheet_index=0; sheet_index < gcmO->ice_regridders().index.size(); ++sheet_index) {

        // Update from ice-only coverage
        {GetSheetElevO sheet(get_sheet_elevO(
            gcmO, paramsO_rawA, sheet_index, emI_ices[sheet_index]));

            for (size_t iO_d=0; iO_d < sheet.dimO.dense_extent(); ++iO_d) {
                auto const iO_s = sheet.dimO.to_sparse(iO_d);
                da_zicetopO(iO_s) += sheet.elev_areaO(iO_d) * gcmO->agridA.native_area(gcmO->agridA.dim.to_dense(iO_s)); //sheet.wO(iO_d);
            }
        }

        {GetSheetElevO sheet(get_sheet_elevO(
            gcmO, paramsO_correctA, sheet_index, emI_ices[sheet_index]));

            for (size_t iO_d=0; iO_d < sheet.dimO.dense_extent(); ++iO_d) {
                auto const iO_s = sheet.dimO.to_sparse(iO_d);
                da_giceO(iO_s) += sheet.wO(iO_d);
            }
        }

        // Update from ice+land coverage
        {GetSheetElevO sheet(get_sheet_elevO(
            gcmO, paramsO_rawA, sheet_index, emI_lands[sheet_index]));

            for (size_t iO_d=0; iO_d < sheet.dimO.dense_extent(); ++iO_d) {
                auto iO_s = sheet.dimO.to_sparse(iO_d);
                da_zatmoO(iO_s) += sheet.elev_areaO(iO_d) * gcmO->agridA.native_area(gcmO->agridA.dim.to_dense(iO_s)); //sheet.wO(iO_d);
            }
        }

        {GetSheetElevO sheet(get_sheet_elevO(
            gcmO, paramsO_correctA, sheet_index, emI_lands[sheet_index]));

            for (size_t iO_d=0; iO_d < sheet.dimO.dense_extent(); ++iO_d) {
                auto iO_s = sheet.dimO.to_sparse(iO_d);
                da_contO(iO_s) += sheet.wO(iO_d);
            }
        }
    }


    // ---------------- Update stuff based on additional land and ice area
    // Update foceanOp, and create new land cells in foceanOm as appropriate.
    for (int iO=0; iO<nO; ++iO) {    // sparse indexing
        if (da_contO(iO) == 0.) continue;

        // Adjust foceanOp, etc. based on how much land we're adding
        double const by_areaO = 1. / gcmO->agridA.native_area(gcmO->agridA.dim.to_dense(iO));

        double  const fgiceOp0 = fgiceOp(iO);
        double const diff_fgiceOp = da_giceO(iO) * by_areaO;
        fgiceOp(iO) = fgiceOp(iO) + diff_fgiceOp;
        foceanOp(iO) = foceanOp(iO) - da_contO(iO)*by_areaO;
        zatmoOp(iO) += da_zatmoO(iO) * by_areaO;

        if (fgiceOp(iO) != 0) {
            double const nzicetopO = (zicetopO(iO)*fgiceOp0 + da_zicetopO(iO)*by_areaO * diff_fgiceOp) / fgiceOp(iO);
            zicetopO(iO) = nzicetopO;
        }

        // When we add more land, some cells that used to be ocean
        // might now become land.
        if ((foceanOp(iO) < 0.5) && (foceanOm(iO) == 1.0)) {
            // Block repeated below
            double const fact = 1. / (1. - foceanOp(iO));
            foceanOm(iO) = 0.0;
            fgiceOm(iO) = fgiceOp(iO) * fact;
            fgrndOm(iO) = 1.0 - fgiceOm(iO) - flakeOm(iO);
            zatmoOm(iO) = zatmoOp(iO) * fact;
        }
    }

    // Remove single-cell oceans in foceanOm (turn more stuff to land)
    long const IM = gcmO->agridA.indexing[0].extent;
    long const JM = gcmO->agridA.indexing[1].extent;
    blitz::TinyVector<int,2> shape2(JM, IM);   // C-style indexing

    for (int j=0; j<JM; ++j) {
    for (int i=0; i<IM; ++i) {

        // Avoid edges, where indexing is more complex (and we don't need to correct anyway)
        if ((i==0) || (i==IM-1) || (j==0) || (j==JM-1)) continue;

        if (foceanOm2(j-1,i) == 0.
            && foceanOm2(j+1,i) == 0.
            && foceanOm2(j,i-1) == 0.
            && foceanOm2(j,i+1) == 0.
            && foceanOm2(j,i) == 1. && foceanOp2(j,i) != 1.)
        {
            long const iO(gcmO->agridA.indexing.tuple_to_index(std::array<long,2>{i,j}));    // Always use alphabetical order

            // Repeated block from above
            double const denom = 1. - foceanOp(iO);
            double const fact = 1. / denom;
            foceanOm(iO) = 0.0;
            fgiceOm(iO) = fgiceOp(iO) * fact;
            fgrndOm(iO) = 1.0 - fgiceOm(iO) - flakeOm(iO);
            zatmoOm(iO) = zatmoOp(iO) * fact;
        }
    }}

    // Run sanity checks on output
    sanity_nonan("foceanOp2", foceanOp2, errors);
    sanity_nonan("fgiceOp2", fgiceOp2, errors);
    sanity_nonan("zatmoOp2", zatmoOp2, errors);
    sanity_nonan("foceanOm2", foceanOm2, errors);
    sanity_nonan("flakeOm2", flakeOm2, errors);
    sanity_nonan("fgrndOm2", fgrndOm2, errors);
    sanity_nonan("fgiceOm2", fgiceOm2, errors);
    sanity_nonan("zatmoOm2", zatmoOm2, errors);
    sanity_nonan("zicetopO2", zicetopO2, errors);
    sanity_check_land_fractions(foceanOm2, flakeOm2, fgrndOm2, fgiceOm2, errors);
}

/** Creates a new indexingHC (1-D Atm. grid indexing plus nhc)
indexing scheme from and old one, changing nhc along the way. */
Indexing indexingHC_change_nhc(Indexing const &indexingHC0, int nhc)
{
    return Indexing(
        std::vector<IndexingData>{
            indexingHC0[0],
            IndexingData(indexingHC0[1].name, indexingHC0[1].base, nhc)
        },
        std::vector<int>(indexingHC0.indices()));
}


EOpvAOpResult compute_EOpvAOp_merged(  // (generates in dense indexing)
SparseSetT &dimAOp,    // dimAOp is appended; dimEOp is returned as part of return variable.
ibmisc::ZArray<int,double,2> const &EOpvAOp_base,    // from linear::Weighted_Compressed; UNSCALED
RegridParams paramsO,
GCMRegridder_Standard const *gcmO,     // A bunch of local ice sheets
double const eq_rad,    // Radius of the earth
std::vector<blitz::Array<double,1>> const &emIs,
bool use_global_ice,
bool use_local_ice,
std::vector<double> const &hcdefs_base, // [nhc]  Elev class definitions for base ice
Indexing const &indexingHC_base,
bool squash_ecs,    // Should ECs be merged if they are the same elevation?
std::vector<std::string> &errors)
{
    EOpvAOpResult ret;    // return variable

    // ======================= Create a merged EOpvAOp of base ice and ice sheets
    // (and then call through to _compute_AAmvEAm)

    // Accumulator for merged EOpvAOp (unscaled)
    MakeDenseEigenT EOpvAOp_m(
        {SparsifyTransform::ADD_DENSE},   // convert sparse to dense indexing
        {&ret.dimEOp, &dimAOp}, '.');
    auto EOpvAOp_accum(EOpvAOp_m.accum());

    // Merge in local matrices
    paramsO.scale = false;
    std::array<long,2> EOpvAOp_sheet_shape {0,0};
    ret.offsetE = 0;
    if (use_local_ice) {
        for (size_t sheet_index=0; sheet_index < gcmO->ice_regridders().index.size(); ++sheet_index) {
            // Get local EOpvAOp matrix
            std::unique_ptr<RegridMatrices_Dynamic> rmO(
                gcmO->regrid_matrices(sheet_index, emIs[sheet_index], paramsO));
            SparseSetT dimEO_sheet, dimAO_sheet;
            std::unique_ptr<ibmisc::linear::Weighted_Eigen> EOpvAOp_sheet(
                rmO->matrix_d("EvA", {&dimEO_sheet, &dimAO_sheet}, paramsO));
            EOpvAOp_sheet_shape = EOpvAOp_sheet->shape();

            // Merge it in...
            // Stack ice sheet ECs on top of global ECs.
            // NOTE: Assumes EC dimension in indexing has largest stride
            for (auto ii(begin(*EOpvAOp_sheet->M)); ii != end(*EOpvAOp_sheet->M); ++ii) {
                // Separate ice sheet ECs from global ECs
                EOpvAOp_accum.add({
                    dimEO_sheet.to_sparse(ii->index(0)),
                    dimAO_sheet.to_sparse(ii->index(1))},
                    ii->value());
            }
        }

        ret.hcdefs.insert(ret.hcdefs.end(), gcmO->hcdefs().begin(), gcmO->hcdefs().end());
        for (size_t i=0; i<gcmO->hcdefs().size(); ++i)
            ret.underice_hc.push_back(UI_ICEBIN);
    }

    std::array<long,2> EOpvAOp_base_shape {0,0};
    if (use_global_ice) {
        ret.offsetE = gcmO->indexingE.extent();
        // Merge in global matrix (compressed format; sparse indexing)
        EOpvAOp_base_shape = EOpvAOp_base.shape();  // sparse shape of ZArray

        // printf("EOpvAOp_base_shape = [%d %d] (NHC=%g)\n", EOpvAOp_base_shape[0], EOpvAOp_base_shape[1], (double)EOpvAOp_base_shape[0] / (double)EOpvAOp_base_shape[1]);

        // Copy elements to accumulator matrix
        for (auto ii(EOpvAOp_base.generator()); ++ii; ) {
            auto iE(ii->index(0));   // In case iE occurs more than once in M

            // Stack global EC's on top of local EC's.
            EOpvAOp_accum.add({ii->index(0)+ret.offsetE, ii->index(1)}, ii->value());
        }
        ret.hcdefs.insert(ret.hcdefs.end(), hcdefs_base.begin(), hcdefs_base.end());
        for (size_t i=0; i<hcdefs_base.size(); ++i)
            ret.underice_hc.push_back(UI_NOTHING);
    }


    // Set overall size
    ret.dimEOp.set_sparse_extent(ret.offsetE + EOpvAOp_base_shape[0]);
    dimAOp.set_sparse_extent(EOpvAOp_base_shape[1]);

    // Convert to Eigen
    ret.EOpvAOp.reset(new EigenSparseMatrixT(EOpvAOp_m.to_eigen()));
    ret.indexingHC = indexingHC_change_nhc(indexingHC_base, ret.hcdefs.size());

    if (squash_ecs) {
        return squash_ECs({&ret.dimEOp, &dimAOp}, *ret.EOpvAOp, ret.hcdefs, ret.indexingHC);
    } else {
        return ret;
    }
}


/** Merges repeated ECs */
EOpvAOpResult squash_ECs(
std::array<SparseSetT *,2> dims0,    // const
EigenSparseMatrixT const &EOpvAOp0,
std::vector<double> const &hcdefs0, // Elevation of each EC in EOpvAOp0
Indexing const &indexingHC0)
{
    EOpvAOpResult ret;

    ConstUniverseT const_dims0({"dimEOp0", "dimAOp0"}, {dims0[0], dims0[1]});

    // Determine new set of ECs
    std::set<double> hcdefs_set(hcdefs0.begin(), hcdefs0.end());
    ret.hcdefs = std::vector<double>(hcdefs_set.begin(), hcdefs_set.end());  // sorted...
    int const nhc = ret.hcdefs.size();

    // Create a map from new HC value to new HC index
    std::map<double,int> hcdefs_map;
    for (int i=0; i<nhc; ++i) {
        hcdefs_map.insert(std::make_pair(ret.hcdefs[i], i));
        ret.underice_hc.push_back(UI_NOTHING);
    }

    // Define to_new[] such that: ihc_new == to_new[ihc_old]
    std::vector<int> to_new;
    for (size_t i=0; i<hcdefs0.size(); ++i)
        to_new.push_back(hcdefs_map[hcdefs0[i]]);

    // Create indexing for for new matrix
    ret.indexingHC = indexingHC_change_nhc(indexingHC0, nhc);


    // Re-do the sparsematrix with new ECs
    MakeDenseEigenT EOpvAOp_m(
        {SparsifyTransform::ADD_DENSE, SparsifyTransform::KEEP_DENSE},   // convert sparse to dense indexing
        std::array<SparseSetT *,2>{&ret.dimEOp, dims0[1]}, '.');
    auto EOpvAOp_accum(EOpvAOp_m.accum());

    std::array<int,2> iTuple;
        int &iAO0(iTuple[0]);
        int &ihc0(iTuple[1]);
    for (auto ii(begin(EOpvAOp0)); ii != end(EOpvAOp0); ++ii) {
        // Convert iEO0 (in old elevation class space) to iE1 (in new)
        int const iEO0 = dims0[0]->to_sparse(ii->index(0));
        indexingHC0.index_to_tuple(&iTuple[0], iEO0);
        int const ihc1 = to_new[ihc0];
        int const iE1 = ret.indexingHC.tuple_to_index(std::array<long,2>{iAO0,ihc1});

        // Separate ice sheet ECs from global ECs
        EOpvAOp_accum.add({iE1, ii->index(1)}, ii->value());
    }
    auto const nA = dims0[1]->sparse_extent();
    ret.dimEOp.set_sparse_extent(nA * nhc);

    ret.EOpvAOp.reset(new EigenSparseMatrixT(EOpvAOp_m.to_eigen()));


    return ret;
}

}}    // namespace
