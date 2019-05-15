#include <cmath>
#include <ibmisc/fortranio.hpp>
#include <ibmisc/ncbulk.hpp>
#include <icebin/error.hpp>
#include <icebin/modele/make_topoo.hpp>
#include <icebin/modele/grids.hpp>

using namespace blitz;
using namespace ibmisc;
using namespace spsparse;

namespace icebin {
namespace modele {


class GetSheetElevO {
    TmpAlloc _tmp;
public:
    SparseSetT dimO;
    blitz::Array<double,1> wO;        // Ice (or land) covered area of each gridcell [m^2]
    blitz::Array<double,1> elevO;     // Elevation [m]
};

GetSheetElevO get_sheet_elevO(
GCMRegridder *gcmO,
RegridParams const &paramsA,
int sheet_index,
ElevMask<1> &emI,
bool include_ice,
bool include_bedrock)
{
    GetSheetElevO ret;

    size_t const nI = gcmO->nI(sheet_index);
    size_t const nO = gcmO->nA();

    // Construct an elevmaskI for ice only
    // (==elevI on ice, nan on ocean or bare land)
    auto &emI(elevmasks[sheet_index]);        // Original mask
    blitz::Array<double,1> elevmaskI(nI);     // Mask in format needed by IceBin
    for (int iI=0; iI<nI; ++iI) {
        switch(emI.mask(iI)) {
            case IceMask::GROUNDED_ICE :
            case IceMask::FLOATING_ICE :
                elevmaskI(iI) = (include_ice ? emI.elev(iI) : nan);
            break;
            case IceMask::ICE_FREE_OCEAN :
            case IceMask::UNKNOWN :
                elevmaskI(iI) = nan;
            break;
            case IceMask::ICE_FREE_BEDROCK :
                elevmask(iI) = (include_bedrock ? emI.elev(iI) : nan);
            break;
        }

    }

    // Obtain the OvI matrix
    SparseSet dimI;
    std::unique_ptr<RegridMatrices_Dynamic> rmO(
        gcmO->regrid_matrices(sheet_index, elevmaskI));
    RegridParams paramsO(paramsA);
        paramsO.scale = false;
        paramsO.correctA = true;
   auto OvI(rmO->matrix_d("AvI", {&ret.dimO, &dimI}, paramsO));

    // Construct elevI in dense indexing space
    ElevMask<1> const &emI(elevmasks[sheet_index]);
    blitz::Array<double,1> elevI(dimI.size());
    for (size_t iI_d=0; iI_d < dimIp.size(); ++iI_d) {
        auto iI_s = dimI.to_sparse(iI_d);
        elevI(iO_d) = emI.elev(iI_s);
    }

    ret.elevO.reference(OvI_ice->apply(elevI_d, nan, true, ret._tmp));    // dense indexing
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
    std::vector<ElevMask<1>> const &elevmasks,    // elevation and cover types for each ice sheet
    double const eq_rad)    // Radius of the earth
{

    auto foceanOp(reshape1(foceanOp2));
    auto fgicepO(reshape1(fgiceOp2));
    auto zatmoOp(reshape1(zatmoOp2));
    auto foceanOm(reshape1(foceanOm2));
    auto fgrndOm(reshape1(fgrndOm2));
    auto fgiceOm(reshape1(fgiceOm2));
    auto zatmoOm(reshape1(zatmoOm2));
    auto zicetopO(reshape1(zicetopO2));

    auto nO = gcmO->nA();

    // -------------- Compute additional area of land and ice due to local ice sheets
    for (size_t sheet_index=0; sheet_index < gcmO->ice_regridders().index.size(); ++sheet_index) {

        // Update from ice-only coverage
        blitz::Array<double,1> da_giceO(nO);     // increase in area of landice [m^2]
        da_giceO = 0.;
        blitz::Array<double,1> da_zicetopO(nO);     // increase in area of landice * elevation [m^3]
        da_zicetopO = 0.;
        {GetSheetElevO sheet(get_sheet_elevO(
            gcmO, paramsA, sheet_index, emI, true, false));

            for (size_t iO_d=0; iO_d < sheet.dimO.size(); ++iO_d) {
                auto const iO_s = sheet.dimO.to_sparse(iO_d);
                da_giceO(iO_s) += sheet.wO(iO_d);
                da_zicetopO(iO_s) += sheet.elevO(iO_d) * sheet.wO(iO_d);
            }
        }

        // Update from ice+land coverage
        blitz::Array<double,1> da_contO(nO);    // increase in area of continent [m^2]
        da_contO = 0.;
        blitz::Array<double,1> da_zatmoO(nO);    // increase in area of continent * elevation [m^3]
        da_zatmoO = 0.;
        {GetSheetElevO sheet(get_sheet_elevO(
            gcmO, paramsA, sheet_index, emI, true, true));

            for (size_t iO_d=0; iO_d < sheet.dimO.size(); ++iO_d) {
                auto iO_s = dimO.to_sparse(iO_d);
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

    auto &foceanOm2(reshape<double,1,2>(foceanOm, shape2));
    for (int j=0; j<JM; ++j) {
    for (int i=0; i<IM; ++i) {

        // Avoid edges, where indexing is more complex (and we don't need to correct anyway)
        if ((i==0) || (i==im-1) || (j==0) || (j==jm-1)) continue;

        if (foceanOm2(j-1,i) == 0.
            && foceanOm2(j+1,i) == 0.
            && foceanOm2(j,i-1) == 0.
            && foceanOm2(j,i+1) == 0.
            && foceanOm2(j,i) == 1.)
        {
            long const iO(gcmO->agridA.indexing.tuple_to_index({j,i}));

            // Repeated block from above
            double const fact = 1. / (1. - foceanOp(iO));
            foceanOm(iO) = 0.0;
            fgrndOm(iO) = 1.0;
            fgiceOm(iO) = fgiceOp(iO) * fact;
            zatmoOm(iO) = zatmoOp(iO) * fact;
        }
    }
}




void merge_topoA(
GCMRegridder_ModelE *gcmA,   // Gets updated with new fcoeanOp, foceanOm
// AAmvEAM is either read from output of global_ec (for just global ice);
// or it's the output of compute_AAmvEAm_merged (for merged global+local ice)
linear::WeightedEigen const &AAmvEAm,

blitz::Array<double,2> const &foceanOm2,     // Rounded FOCEAN
blitz::Array<double,2> const &flakeOm2,
blitz::Array<double,2> const &fgrndOm2,
blitz::Array<double,2> const &fgiceOm2,
blitz::Array<double,2> const &zatmoOm2,
blitz::Array<double,2> const &zlakeOm2,
blitz::Array<double,2> const &zicetopOm2,
//SparseSetT const &dimO,    // Tells us which grid cells in O were changed.
blitz::Array<double,2> &foceanA2,    // Rounded FOCEAN
blitz::Array<double,2> &flakeA2,
blitz::Array<double,2> &fgrndA2,
blitz::Array<double,2> &fgiceA2,
blitz::Array<double,2> &zatmoA2,
blitz::Array<double,2> &zlakeA2,
blitz::Array<double,2> &zicetopA2)
//
EigenSparseMatrixT const &EOpvAOp,
RegridParams const &paramsA,
SparseSetT const &dimEOp,
SparseSetT const &dimAOp,
double const eq_rad,
//
blitz::Array<double,3> &fhc3,
blitz::Array<double,3> &elevE3,
blitz::Array<uint16_t,3> &underice3
{
    auto foceanOm(reshape1(foceanOm2));
    auto flakeO(reshape1(flakeO2));
    auto fgiceOm(reshape1(fgiceOm2));
    auto zatmoOm(reshape1(zatmoOm2));
    auto zicetopO(reshape1(zicetopO2));

    auto foceanA(reshape1(foceanA2));
    auto flakeA(reshape1(flakeA2));
    auto fgrndA(reshape1(fgrndA2));
    auto fgiceA(reshape1(fgiceA2));
    auto zatmoA(reshape1(zatmoA2));

    gcmA->foceanOp = reshape1(foceanOp);    // COPY
    gcmA->foceanOm = reshape1(foceanOm);    // COPY


    HntrSpec &hntrO(dynamic_cast<GridSpec_LonLat *>(&*gcmA->gcmO->agridA.spec)->hntr);
    HntrSpec &hntrA(dynamic_cast<GridSpec_LonLat *>(&*gcmA->agridA.spec)->hntr);
    Hntr hntr_AvO(17.17, meta.hspecA, hspecO);


    typedef std::array<blitz::Array<double,2> &,3> tup;
    std::vector<oatup> tups {
        {foceanOm2, foceanA2, WTO},
        {flakeOm2, flakeA2, WTO},
        {fgrndOm2, fgrndA2, WTO},
        {fgiceOm2, fgiceA2, WTO},
        {zatmoOm2, zatmoA2, WTO},
        {zlakeOm2, zlakeA2, WTO},
        {zicetopOm2, zicetopA2, fgiceO}};

    // Regrid TOPOO variables and save to TOPOA
    blitz::Array<double, 2> WTO(const_array(blitz::shape(hspecO.jm,hspecO.im), 1.0));
    for (auto &t : tups) {
        auto &xxOm2(t[0]);
        auto &xxA2(t[1]);
        auto &xxWTO(t[2]);

        // Regrid to A grid
        hntr_AvO.regrid(xxWTO, xxOm2, xxA2);
    }

    // ================= Create fhc, elevE and underice
    int const nhc_icebin = gcmA->nhc();  // = gcmA->indexingE[2].extent
    int const nhc_gcm = 1 + nhc_icebin;
    blitz::TinyVector<int,3> shapeE(nhc_gcm, gcmA->indexingE[1].extent, gcmA->indexingE[0].extent);
    blitz::TinyVector<int,2> shapeE2(shapeE(0), shapeE(1)*shapeE(2));

    // Initialize
    fhc = 0;
    elevE = NaN;
    underice = 0;    // Elevation grid cell unused

    auto all(blitz::Range::all());

    // ------------ Segment 0: Elevation Classes
    printf("Segment 2: Elevation Classes\n");
    auto fhcE2(reshape<double,3,2>(fhc3, shapeE2));
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
        EOpvAOp, dimEOp, dimAOp, wAOp);

    for (auto ii=begin(AAmvEAm); ii != end(AAmvEAm); ++ii) {
        auto const iA_d = ii->index(0);
        auto const iA = dimAAm.to_sparse(iA_d);
        auto const iE_d = ii->index(1);
        auto const iE = dimEAm.to_sparse(iE_d);

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

        fhcE2(ihc,iA) += ii->value() / AAmvEAm->wM(iA);   // AAMvEAM is not scaled
        undericeE2(ihc,iA) = gcmA->underice(ihc);UI_NOTHING;    // No IceBin coupling here
        elevE2(iA) = gcmA->hcdefs[ihc];
    }

    // ------------ Segment 1: land part of sealand
    printf("Segment 1: seaLAND\n");
    ec_base = nhc_icebin;    
    for (int j=0; j<meta.hspecA.jm; ++j) {
    for (int i=0; i<meta.hspecA.im; ++i) {
        if (fgiceA(j,i) > 0) {
            fhc(ec_base, j,i) = 1e-30;
            elevE(ec_base, j,i) = zicetopA(j,i);
            underice(ec_base, j,i) = UI_NOTHING;
        }
    }}


}}
