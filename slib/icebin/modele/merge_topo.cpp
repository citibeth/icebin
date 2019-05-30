#include <icebin/modele/merge_topo.hpp>
#include <icebin/eigen_types.hpp>
#include <ibmisc/linear/compressed.hpp>

using namespace blitz;
using namespace ibmisc;
using namespace spsparse;

static double const NaN = std::numeric_limits<double>::quiet_NaN();

namespace icebin {
namespace modele {

static blitz::Array<double,1> get_elevmaskI(ElevMask<1> const &emI, bool include_ice, bool include_bedrock)
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

class GetSheetElevO {
public:
    TmpAlloc _tmp;
    SparseSetT dimO;
    blitz::Array<double,1> wO;        // Ice (or land) covered area of each gridcell [m^2]
    blitz::Array<double,1> elevO;     // Elevation [m]
};

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
GCMRegridder const *gcmO,     // A bunch of local ice sheets
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

}}    // namespace
