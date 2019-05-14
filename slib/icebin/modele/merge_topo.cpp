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
    blitz::Array<double,2> const &areaO,    // Area of each ocean gridcell
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
blitz::Array<double,2> const &flakeO2,
blitz::Array<double,2> const &fgrndO2,
blitz::Array<double,2> const &fgiceO2,
blitz::Array<double,2> const &zatmoO2,
//SparseSetT const &dimO,    // Tells us which grid cells in O were changed.
blitz::Array<double,2> const &foceanA2,    // Rounded FOCEAN
blitz::Array<double,2> const &flakeA2,
blitz::Array<double,2> const &fgrndA2,
blitz::Array<double,2> const &fgiceA2,
blitz::Array<double,2> const &zatmoA2)
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

    // =====================================================
    // Regrid TOPO to Atmosphere grid
    HntrSpec const &hntrA(cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr);
    HntrSpec const &hntrO(cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr);
    Hntr hntrAvO(17.17, hntrA, hntrO);

    TupleListT<2> AvO_tp;
    hntrAvO.scaled_regrid_matrix(spsparse::accum::ref(AvO_tp));
    EigenSparseMatrixT AvO_e(hntrA.size(), hntrO.size());
    AvO_e.setFromTriplets(AvO_tp.begin(), AvO_tp.end());

    fgiceA = 0;

    map_eigen_colvector(foceanA) = AvO_e * map_eigen_colvector(foceanOm);
    map_eigen_colvector(flakeA) = AvO_e * map_eigen_colvector(flakeO);
    map_eigen_colvector(fgrndA) = AvO_e * map_eigen_colvector(fgrndO);
    map_eigen_colvector(fgiceA) = AvO_e * map_eigen_colvector(fgiceO);
    map_eigen_colvector(zatmoA) = AvO_e * map_eigen_colvector(zatmoO);
    map_eigen_colvector(elevA) = AvO_e * map_eigen_colvector(elevO);


    TupleListT<2> AvE_global_tp;
    blitz::Array<double,1> elevE_global(nE);


}




void merge_topo_O_A()
{
}

