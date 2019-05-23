#include <ibmisc/linear/compressed.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/modele/GCMRegridder_ModelE.hpp>
#include <icebin/modele/hntr.hpp>

using namespace blitz;
using namespace ibmisc;
using namespace spsparse;

static double const NaN = std::numeric_limits<double>::quiet_NaN();

namespace icebin {
namespace modele {


class GetSheetElevO {
public:
    TmpAlloc _tmp;
    SparseSetT dimO;
    blitz::Array<double,1> wO;        // Ice (or land) covered area of each gridcell [m^2]
    blitz::Array<double,1> elevO;     // Elevation [m]
};

blitz::Array<double,1> get_elevmaskI(ElevMask<1> const &emI, bool include_ice, bool include_bedrock)
{
    auto nI(emI.elev.extent(0));

    blitz::Array<double,1> elevmaskI(nI);     // Mask in format needed by IceBin
    for (int iI=0; iI<nI; ++iI) {
        switch(emI.mask(iI)) {
            case IceMask::GROUNDED_ICE :
            case IceMask::FLOATING_ICE :
                elevmaskI(iI) = (include_ice ? emI.elev(iI) : NaN);
            break;
            case IceMask::ICE_FREE_OCEAN :
            case IceMask::UNKNOWN :
                elevmaskI(iI) = NaN;
            break;
            case IceMask::ICE_FREE_BEDROCK :
                elevmaskI(iI) = (include_bedrock ? emI.elev(iI) : NaN);
            break;
        }
    }
    return elevmaskI;
}

GetSheetElevO get_sheet_elevO(
GCMRegridder *gcmO,
RegridParams const &paramsA,
int sheet_index,
std::vector<ElevMask<1>> const &elevmasks,    // elevation and cover types for each ice sheet
bool include_ice,
bool include_bedrock)
{
    GetSheetElevO ret;

    size_t const nI = gcmO->nI(sheet_index);
    size_t const nO = gcmO->nA();

    // Construct an elevmaskI for ice only
    // (==elevI on ice, NaN on ocean or bare land)
    ElevMask<1> const &emI(elevmasks[sheet_index]);        // Original mask
    blitz::Array<double,1> elevmaskI(get_elevmaskI(emI, include_ice, include_bedrock));

    // Obtain the OvI matrix
    SparseSetT dimI;
    std::unique_ptr<RegridMatrices_Dynamic> rmO(
        gcmO->regrid_matrices(sheet_index, elevmaskI));
    RegridParams paramsO(paramsA);
        paramsO.scale = false;
        paramsO.correctA = true;
   auto OvI(rmO->matrix_d("AvI", {&ret.dimO, &dimI}, paramsO));

    // Construct elevI in dense indexing space
    blitz::Array<double,1> elevI(dimI.dense_extent());
    for (size_t iI_d=0; iI_d < dimI.dense_extent(); ++iI_d) {
        auto iI_s = dimI.to_sparse(iI_d);
        elevI(iI_d) = emI.elev(iI_s);
    }

    ret.elevO.reference(OvI->apply(elevI, NaN, true, ret._tmp));    // dense indexing
    ret.wO.reference(OvI->wM);

    return ret;
}

void merge_topoO(
    // ------ TOPOO arrays, originally with just global ice
    // Ice model viewpoint (fractional ocean cells)
    blitz::Array<double,2> &foceanOp2,    // Fractional FOCEAN
    blitz::Array<double,2> &fgiceOp2,
    blitz::Array<double,2> &zatmoOp2,
    // ModelE viewpoint (rounded ocean cells)
    blitz::Array<double,2> &foceanOm2,     // Rounded FOCEAN
    blitz::Array<double,2> &fgrndOm2,
    blitz::Array<double,2> &fgiceOm2,
    blitz::Array<double,2> &zatmoOm2,
    // Not affected by Om; as long as top of ice is maintained even for ocean-rounded cells.
    blitz::Array<double,2> &zicetopO2,
    // ------ Local ice to merge in...
    GCMRegridder *gcmO,
    RegridParams const &paramsA,
    std::vector<ElevMask<1>> const &elevmasks,    // elevation and cover types for each ice sheet
    double const eq_rad)    // Radius of the earth
{

    auto foceanOp(reshape1(foceanOp2));
    auto fgiceOp(reshape1(fgiceOp2));
    auto zatmoOp(reshape1(zatmoOp2));
    auto foceanOm(reshape1(foceanOm2));
    auto fgrndOm(reshape1(fgrndOm2));
    auto fgiceOm(reshape1(fgiceOm2));
    auto zatmoOm(reshape1(zatmoOm2));
    auto zicetopO(reshape1(zicetopO2));

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

    // -------------- Compute additional area of land and ice due to local ice sheets
    for (size_t sheet_index=0; sheet_index < gcmO->ice_regridders().index.size(); ++sheet_index) {

        // Update from ice-only coverage
        {GetSheetElevO sheet(get_sheet_elevO(
            gcmO, paramsA, sheet_index, elevmasks, true, false));

            for (size_t iO_d=0; iO_d < sheet.dimO.dense_extent(); ++iO_d) {
                auto const iO_s = sheet.dimO.to_sparse(iO_d);
                da_giceO(iO_s) += sheet.wO(iO_d);
                da_zicetopO(iO_s) += sheet.elevO(iO_d) * sheet.wO(iO_d);
            }
        }

        // Update from ice+land coverage
        {GetSheetElevO sheet(get_sheet_elevO(
            gcmO, paramsA, sheet_index, elevmasks, true, true));

            for (size_t iO_d=0; iO_d < sheet.dimO.dense_extent(); ++iO_d) {
                auto iO_s = sheet.dimO.to_sparse(iO_d);
                da_contO(iO_s) += sheet.wO(iO_d);
                da_zatmoO(iO_s) += sheet.elevO(iO_d) * sheet.wO(iO_d);
            }
        }
    }


    // ---------------- Update stuff based on additional land and ice area
    // Update foceanOp, and create new land cells in foceanOm as appropriate.
    for (int iO=0; iO<nO; ++iO) {
        if (da_contO(iO) == 0.) continue;

        // Adjust foceanOp, etc. based on how much land we're adding
        double const by_areaO = 1. / gcmO->agridA.native_area(gcmO->agridA.dim.to_dense(iO));

        double  const fgiceOp0 = fgiceOp(iO);
        fgiceOp(iO) = std::min(1.0, fgiceOp(iO) + da_giceO(iO) * by_areaO);
        foceanOp(iO) = std::max(0.0, foceanOp(iO) - da_contO(iO)*by_areaO);
        zatmoOp(iO) += da_zatmoO(iO) * by_areaO;

        zicetopO(iO) = (zicetopO(iO)*fgiceOp0 + da_zicetopO(iO)*by_areaO)
            / (fgiceOp0 + da_giceO(iO)*by_areaO);

        // When we add more land, some cells that used to be ocean
        // might now become land.
        if ((foceanOp(iO) < 0.5) && (foceanOm(iO) == 1.0)) {
            // Block repeated below
            double const fact = 1. / (1. - foceanOp(iO));
            foceanOm(iO) = 0.0;
            fgrndOm(iO) = 1.0;
            fgiceOm(iO) = fgiceOp(iO) * fact;
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
            && foceanOm2(j,i) == 1.)
        {
            long const iO(gcmO->agridA.indexing.tuple_to_index(std::array<long,2>{j,i}));

            // Repeated block from above
            double const fact = 1. / (1. - foceanOp(iO));
            foceanOm(iO) = 0.0;
            fgrndOm(iO) = 1.0;
            fgiceOm(iO) = fgiceOp(iO) * fact;
            zatmoOm(iO) = zatmoOp(iO) * fact;
        }
    }}
}



EigenSparseMatrixT compute_EOpvAOp_merged(  // (generates in dense indexing)
    std::array<SparseSetT *,2> dims,
    std::string const &global_ecO,    // File written by global_ec
    RegridParams const &paramsA,
    GCMRegridder const *gcmO,
    double const eq_rad,    // Radius of the earth
    std::vector<ElevMask<1>> const elevmasks,    // elevation and cover types for each ice sheet
    bool use_global_ice,
    bool use_local_ice,
    bool include_bedrock)    // true if non-ice covered areas of land should also be included
{
    // ======================= Create a merged EOpvAOp of base ice and ice sheets
    // (and then call through to _compute_AAmvEAm)

    // Accumulator for merged EOpvAOp (unscaled)
    MakeDenseEigenT EOpvAOp_m(
        {SparsifyTransform::ADD_DENSE},   // convert sparse to dense indexing
        dims, '.');

    // Regrid params for ice sheet regrid matrices
    RegridParams paramsO(paramsA);
        paramsO.scale = false;
        paramsO.correctA = true;

    // Merge in local matrices
    std::array<long,2> EOpvAOp_sheet_shape {0,0};
    if (use_local_ice) {
        for (size_t sheet_index=0; sheet_index < gcmO->ice_regridders().index.size(); ++sheet_index) {
            ElevMask<1> const &emI(elevmasks[sheet_index]);
            int const nI(emI.elev.extent(0));

            // Construct an elevmaskI for ice sheet, =NaN off ice sheet
            blitz::Array<double,1> elevmaskI(get_elevmaskI(emI, true, include_bedrock));

            // Get local EOpvAOp matrix
            std::unique_ptr<RegridMatrices_Dynamic> rmO(
                gcmO->regrid_matrices(sheet_index, elevmaskI, paramsO));
            SparseSetT dimEO_sheet, dimAO_sheet;
            std::unique_ptr<ibmisc::linear::Weighted_Eigen> EOpvAOp_sheet(
                rmO->matrix_d("EvA", {&dimEO_sheet, &dimAO_sheet}, paramsO));
            EOpvAOp_sheet_shape = EOpvAOp_sheet->shape();

            // Merge it in...
            // Stack ice sheet ECs on top of global ECs.
            // NOTE: Assumes EC dimension in indexing has largest stride
            for (auto ii(begin(*EOpvAOp_sheet->M)); ii != end(*EOpvAOp_sheet->M); ++ii) {
                // Separate ice sheet ECs from global ECs
                EOpvAOp_m.M.add({
                    dimEO_sheet.to_sparse(ii->index(0)),
                    dimAO_sheet.to_sparse(ii->index(1))},
                    ii->value());
            }
        }
    }

    std::array<long,2> EOpvAOp_base_shape {0,0};
    long offsetE = gcmO->indexingE.extent();
    if (use_global_ice) {
        // Merge in global matrix (compressed format; sparse indexing) If
        // an EC is already used by an ice sheet, remove the ice in that
        // EC from the global matrix.  This "cheap trick" will result in a
        // SMALL amount of missing ice.  But it's simpler than creating a
        // new set of EC's for global vs. ice sheet ice
        {NcIO ncio(global_ecO, 'r');
            // Load from the file
            linear::Weighted_Compressed EOpvAOp_base;
            EOpvAOp_base.ncio(ncio, "EvA");  // (A in this file) == O
            std::array<long,2> EOpvAOp_base_shape(EOpvAOp_base.shape());

            // Check that global matrix has same number of ECs as local
            // (and hopefully at same elevations too)
            if (gcmO->nhc() != EOpvAOp_base_shape[0] / EOpvAOp_base_shape[1])
                (*icebin_error)(-1, "NHC mismatch between global and local ice");

            // Copy elements to accumulator matrix
            for (auto ii(EOpvAOp_base.M.generator()); ++ii; ) {
                auto iE(ii->index(0));   // In case iE occurs more than once in M

                // Stack global EC's on top of local EC's.
                EOpvAOp_m.M.add({ii->index(0)+offsetE, ii->index(1)}, ii->value());
            }
        }
    }
        
    if ((EOpvAOp_sheet_shape[0] != EOpvAOp_base_shape[0]) ||
        (EOpvAOp_sheet_shape[1] != EOpvAOp_base_shape[1])) {

        (*icebin_error)(-1, "Shape of on-disk and local EOpvAOp differs (%ld,%ld) vs (%ld,%ld)!",
            EOpvAOp_base_shape[0], EOpvAOp_base_shape[1],
            EOpvAOp_sheet_shape[0], EOpvAOp_sheet_shape[1]);
    }

    // Set overall size
    dims[0]->set_sparse_extent(offsetE + EOpvAOp_base_shape[0]);
    dims[1]->set_sparse_extent(EOpvAOp_base_shape[1]);

    // Compute EOpvAOp and wAOp
    EigenSparseMatrixT EOpvAOp(EOpvAOp_m.to_eigen());
}



std::unique_ptr<linear::Weighted_Eigen> _compute_AAmvEAm(
    std::array<SparseSetT *,2> dims,
    RegridParams const &paramsA,
    double const eq_rad,    // Radius of the earth

    // Things obtained from gcmA
    unsigned long const nA,     // gcmA->nA()
    unsigned long const nE,    // gcmA->nE()
    unsigned int const nhc,    // gcmA->nhc()
    HntrSpec const &hntrO,        // cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr
    HntrSpec const &hntrA,        // cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr
    IndexSet const indexingHCO,    // gcmA->gcmO->indexingHC
    IndexSet const indexingHCA,    // gcmA->indexingHC
    blitz::Array<double,1> const &foceanAOp,    // gcmA->foceanAOp
    blitz::Array<double,1> const &foceanAOm,    // gcmA->foceanAOm

    // Sub-parts of the computation, pre-computed
    EigenSparseMatrixT const &EOpvAOp,
    SparseSetT &dimEOp,
    SparseSetT &dimAOp,
    blitz::Array<double,1> const &wAOp)
{

    ConstUniverse const_dimAOp({"dimEOp", "dimAOp"}, {&dimEOp, &dimAOp});

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
    if (paramsA.scale) {
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
GCMRegridder_ModelE *gcmA,   // Gets updated with new fcoeanOp, foceanOm
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
EigenSparseMatrixT const &EOpvAOp,
RegridParams const &paramsA,
SparseSetT &dimEOp,    // const
SparseSetT &dimAOp,    // const
double const eq_rad,
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
    ConstUniverse const_dimAOp({"dimEOp", "dimAOp"}, {&dimEOp, &dimAOp});

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

    gcmA->foceanAOp = reshape1(foceanOp2);    // COPY
    gcmA->foceanAOm = reshape1(foceanOm2);    // COPY


    HntrSpec &hspecO(dynamic_cast<GridSpec_LonLat *>(&*gcmA->gcmO->agridA.spec)->hntr);
    HntrSpec &hspecA(dynamic_cast<GridSpec_LonLat *>(&*gcmA->agridA.spec)->hntr);
    Hntr hntr_AvO(17.17, hspecA, hspecO);


    blitz::Array<double, 2> WTO(const_array(blitz::shape(hspecO.jm,hspecO.im), 1.0));
    hntr_AvO.regrid(WTO, foceanO2, foceanA2);
    hntr_AvO.regrid(WTO, flakeO2, flakeA2);
    hntr_AvO.regrid(WTO, fgrndO2, fgrndA2);
    hntr_AvO.regrid(WTO, fgiceO2, fgiceA2);
    hntr_AvO.regrid(WTO, zatmoO2, zatmoA2);
    hntr_AvO.regrid(WTO, zlakeO2, zlakeA2);
    hntr_AvO.regrid(fgiceO, zicetopO2, zicetopA2);

    // ================= Create fhc, elevE and underice
    int const nhc_icebin = gcmA->nhc();  // = gcmA->indexingE[2].extent
    int const nhc_gcm = 1 + nhc_icebin;
    blitz::TinyVector<int,3> shapeE(nhc_gcm, gcmA->indexingE[1].extent, gcmA->indexingE[0].extent);
    blitz::TinyVector<int,2> shapeE2(shapeE(0), shapeE(1)*shapeE(2));

    // Initialize
    fhc3 = 0;
    elevE3 = NaN;
    underice3 = 0;    // Elevation grid cell unused

    auto all(blitz::Range::all());

    // ------------ Segment 0: Elevation Classes
    printf("Segment 2: Elevation Classes\n");
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
        {&dimAAm, &dimEAm}, paramsA, gcmA, eq_rad,
        EOpvAOp, dimEOp, dimAOp, wAOp));

    for (auto ii=begin(*AAmvEAm->M); ii != end(*AAmvEAm->M); ++ii) {
        int const iA_d = ii->index(0);
        int const iA = dimAAm.to_sparse(iA_d);
        int const iE_d = ii->index(1);
        int const iE = dimEAm.to_sparse(iE_d);

        gcmA->indexingHC.index_to_tuple(&iTuple[0], iE);

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
        fhcE2(ihc,iA) += ii->value() / AAmvEAm->wM(iA);    // AAmvEAm is not scaled
        undericeE2(ihc,iA) = gcmA->underice(ihc);UI_NOTHING;    // No IceBin coupling here
        elevE2(iA) = gcmA->hcdefs()[ihc];
    }

    // ------------ Segment 1: land part of sealand
    printf("Segment 1: seaLAND\n");
    int ec_base = nhc_icebin;    
    for (int j=0; j<hspecA.jm; ++j) {
    for (int i=0; i<hspecA.im; ++i) {
        if (fgiceA(j,i) > 0) {
            fhc3(ec_base, j,i) = 1e-30;
            elevE3(ec_base, j,i) = zicetopA2(j,i);
            underice3(ec_base, j,i) = UI_NOTHING;
        }
    }}
}













// -----------------------------------------------------------
std::unique_ptr<GCMRegridder> new_gcmA_standard(
    HntrSpec const &hspecA,
    std::string const &grid_name,
    ParseArgs const &args, blitz::Array<double,2> const &elevmaskI)
{
    ExchangeGrid aexgrid;    // Put our answer in here

    auto const &hspecI(args.hspecI);
    modele::Hntr hntr(17.17, hspecA, hspecI);


    // -------------------------------------------------------------
    printf("---- Computing overlaps\n");

    // Compute overlaps for cells with ice
    SparseSet<long,int> _dimA;    // Only include A grid cells with ice
    SparseSet<long,int> _dimI;    // Only include I grid cells with ice
    hntr.overlap(ExchAccum(aexgrid, reshape1(elevmaskI), _dimA, _dimI), args.eq_rad);

    // -------------------------------------------------------------
    printf("---- Creating gcmA for %s\n", grid_name.c_str());

    // Turn HntrSpec --> GridSpec
    GridSpec_LonLat specA(make_grid_spec(hspecA, false, 1, args.eq_rad));
    GridSpec_LonLat specI(make_grid_spec(hspecI, false, 1, args.eq_rad));

    // Realize A grid for relevant gridcells
    auto agridA(make_abbr_grid(grid_name, specA, std::move(_dimA)));

    // Set up elevation classes    
    std::vector<double> hcdefs;
    for (double elev=args.ec_range[0]; elev <= args.ec_range[1]; elev += args.ec_skip) {
        hcdefs.push_back(elev);
    }

    // Create standard GCMRegridder for A <--> I
    std::unique_ptr<GCMRegridder_Standard> gcmA(new GCMRegridder_Standard);
    gcmA->init(
        std::move(agridA), std::move(hcdefs),
        Indexing({"A", "HC"}, {0,0}, {agridA.dim.sparse_extent(), hcdefs.size()}, {1,0}),
        args.correctA);


    // --------------------------------------------------
    // Create IceRegridder for I and add to gcmA
    auto ice(new_ice_regridder(IceRegridder::Type::L0));
    auto agridI(make_abbr_grid("Ice", specI, std::move(_dimI)));
    ice->init("globalI", gcmA->agridA, nullptr,
        std::move(agridI), std::move(aexgrid),
        InterpStyle::Z_INTERP);    // You can use different InterpStyle if you like.

    gcmA->add_sheet(std::move(ice));

    return std::unique_ptr<GCMRegridder>(gcmA.release());
}

std::unique_ptr<GCMRegridder> new_gcmA_mismatched(
    FileLocator const &files, ParseArgs const &args, blitz::Array<double,2> const &elevmaskI)
{
    auto const &hspecO(args.hspecO);
    auto const &hspecI(args.hspecI);

    auto gcmO(new_gcmA_standard(hspecO, "Ocean", args, elevmaskI));


    // --------------------------------------------------
    printf("---- Creating gcmA\n");

    // Create a mismatched regridder, to mediate between different ice
    // extent of GCM vs. IceBin
    std::unique_ptr<modele::GCMRegridder_ModelE> gcmA(
        new modele::GCMRegridder_ModelE("",
            std::shared_ptr<GCMRegridder>(gcmO.release())));

    HntrSpec const &hspecA(cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr);

    // Load the fractional ocean mask (based purely on ice extent)
    {auto fname(files.locate(args.topoo_fname));

        blitz::Array<double,2> foceanO(hspecO.jm, hspecO.im);    // called FOCEAN in make_topoo
        blitz::Array<double,2> foceanfO(hspecO.jm, hspecO.im);    // called FOCEANF in make_topoo

        printf("---- Reading FOCEAN: %s\n", fname.c_str());
        NcIO ncio(fname, 'r');
        ncio_blitz(ncio, foceanO, "FOCEAN", "double", {});
        ncio_blitz(ncio, foceanfO, "FOCEANF", "double", {});


        gcmA->foceanAOp = reshape1(foceanfO);  // COPY: FOCEANF 
        gcmA->foceanAOm = reshape1(foceanO);   // COPY: FOCEAN
    }

    return std::unique_ptr<GCMRegridder>(gcmA.release());
}


}}    // naespace icebin::modele
