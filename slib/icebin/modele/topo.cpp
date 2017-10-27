#include <spsparse/eigen.hpp>

static double const nan = std::numeric_limits<double>::quiet_NaN();

/**

@param fcontI continent mask, indexed by ice sheet ID;
     equal to 1.0 for continent land (or
     ice), 0.0 for ocean.
   NOTE: This is NOT the same as !isnan(elevI), since elevI is
      defined only for the ice sheet.
*/
class ComputeFOCEAN {
    blitz::Array<double,1> foceanOp;
    blitz::Array<double,1> foceanOm;
}


template<int RANK>
struct ElevMask {
    blitz::Array<double,RANK> elev;
    blitz::Array<char,RANK> mask;

    ElevMask(
        blitz::Array<double,RANK> const &_elev,
        blitz::Array<char,RANK> const &_mask)
    : elev(_elev), mask(_mask) {}
};


/** Adds ice sheet information to an FOCEAN read from Gary's topo files.

@param foceanOp foceanO as read in from TOPO files, with ice sheets
    removed.  Starting out, foceanOp should be 0 or 1 everywhere.
    This function will change foceanOp in areas of ice sheets, setting
    to values in the range [0,1]
@param fcont_elevIs Elevation on each grid in the GCMRegridder.
    Should be nan where there is no land; and either 0, or an actual
    elevation, where there is.
*/
void update_foceanOp(
GCMRegridder *gcmO,
std::map<std::string, ElevMask<1>> const &elevmasks,
blitz::Array<double,1> &foceanOp,    // OUT: 0-based array
blitz::Array<char,1> &changedO)    // OUT
{
//    GCMRegridder *gcmA = &*gcmc->gcm_regridder;


    auto nO(gcmO.nA());

    // --------------------- Compute fcontOp_d (and foceanOp)
    for (auto sheet_name=sheets_index.begin(); sheet_name != sheets_index.end(); ++sheet_name) {

        // Construct an elevmaskI for CONTINENTAL land (not just ice sheet)
        // (==elevI on continent, nan on ocean)
        auto &emI(elevmasks.at(*sheet_name));
        int nI = emI.elev.extent(0);
        blitz::Array<double,1> elevmaskI(nI);
        for (int iI=0; iI<nI; ++iI) {
            auto const m(emI.mask(iI));
            elevmaskI(iI) = (m == MASK_ICE_FREE_OCEAN ? nan : emI.elev(iI));
        }

        // Get OvI for continental cells
        RegridMatrices rm(gcmO.regrid_matrices(*sheet_name, elevmaskI));
        SparseSetT dimO, dimI;
        RegridMatrices::Params paramsO;
            paramsO.scale = false;
            paramsO.correctA = false;
        auto OvI(rmO.matrix('AvI', {&dimO, &dimI}, paramsO);

        // Don't need to set up the mask on I ourselves; this is already
        // built into the OvI matrix.  The mask, taken from PISM, includes
        // all bare land and ice-covered areas.
        // See: pygiss/giss/pism.py   _get_landmask()
        //    (used by write_icebin_in_base.py)
        blitz::Array<double,1> fcontI_d(dimI.dense_extent());
        fcontI_d = 1.0;

        // Compute fcontOp (for this ice sheet only)
        blitz::Array<double,1> fcontOp_d(OvI.apply(fcontI_d, 0., true));    // force_conservation set to true by default, and it probably doesn't matter; but I think it should be false here.

        // Interpolate into foceanOp_s
        for (int iO_d=0; iO_d<fcontOp_d.extent(0); ++iO_d) {
            auto const iO_s = dimO.to_sparse(iO_d);
            foceanOp(iO_s) -= fcontOp_d(iO_d);
            changedO(iO_s) = 1;    // true
        }
    }
}

/** Adds ice sheets to Gary's TOPO file */
void update_fgiceO(
GCMRegridder *gcmO,
std::map<std::string, ElevMask<1>> const &elevmasks,
blitz::Array<double,1> &fgiceO,    // OUT: 0-based array
blitz::Array<char,1> &changedO)    // OUT
{
    auto nO(gcmO.nA());

    // --------------------- Compute fgiceO
    for (auto sheet_name=sheets_index.begin(); sheet_name != sheets_index.end(); ++sheet_name) {

        // Construct an elevmaskI for CONTINENTAL land and ICE SHEET
        // (==elevI on continent, nan on ocean)
        auto &emI(elevmasks.at(*sheet_name));
        int nI = emI.elev.extent(0);
        blitz::Array<double,1> elevmaskI(nI);
        for (int iI=0; iI<nI; ++iI) {
            auto const m(emI.mask(iI));
            elevmask(iI) = (
                m==MASK_GROUNDED_ICE || m==MASK_FLOATING_ICE ?
                emI.elev(iI) : nan);
        }

        // Get OvI for ice cells
        RegridMatrices rm(gcmO.regrid_matrices(*sheet_name, elevmaskI));
        SparseSetT dimO, dimI;
        RegridMatrices::Params paramsO;
            paramsO.scale = false;
            paramsO.correctA = false;
        auto OvI(rmO.matrix('AvI', {&dimO, &dimI}, paramsO);

        // Don't need to set up the mask on I ourselves; this is already
        // built into the OvI matrix.  The mask, taken from PISM, includes
        // all bare land and ice-covered areas.
        // See: pygiss/giss/pism.py   _get_landmask()
        //    (used by write_icebin_in_base.py)
        blitz::Array<double,1> fgiceI_d(dimI.dense_extent());
        fgiceI_d = 1.0;

        // Compute fgiceO (for this ice sheet only)
        blitz::Array<double,1> fgiceO_d(OvI.apply(fgiceI_d, 0., true));    // force_conservation set to true by default, and it probably doesn't matter; but I think it should be false here.

        // Interpolate into foceanOp_s
        for (int iO_d=0; iO_d<fgiceO_d.extent(0); ++iO_d) {
            auto const iO_s = dimO.to_sparse(iO_d);
            fgiceO(iO_s) += fgiceO_d(iO_d);
            changedO(iO_s) = 1;    // true
        }
    }
}


// ======================================================================

/** This needs to be run at least once before matrices can be generated. */
void GCMCoupler_ModelE::update_topo(double time_s, bool first)
{
    GCMRegridder_ModelE const *gcmA = dynamic_cast<GCMRegridder_ModelE *>(&*gcm_regridder);

    std::map<std::string, ElevMask<1>> const &elevmasks;

    for (auto sheet_name=sheets_index.begin(); sheet_name != sheets_index.end(); ++sheet_name) {
        IceCoupler *icec(ice_coupler(*sheet_name));
        elevmasks.insert(std::make_pair(*sheet_name, ElevMask<1>(icec->elevI, icec->maskI)));
    }

    Topos &topoA(modele_inputs);
    update_topo(gcmA, topoO_fname, first,
        topoA, foceanOm0);
}


void update_topo(
    GCMRegridder_ModelE *gcmA,    // Gets updated with new fcoeanOp, foceanOm
    std::string const &topoO_fname,    // Name of Ocean-based TOPO file (aka Gary)
    std::map<std::string, ElevMask<1>> const &elevmasks,
    bool first,    // true if this is the first (initialization) timestep
    HCSegmentData con
    // OUTPUT parameters (variables come from GCMCoupler); must be pre-allocated
    Topos &topoA,
    biltz::Array<double,1> foceanOm0)
{
    if (!first) (*icebin_error)(-1,
        "GCMCoupler_ModelE::update_topo() currently only works for the initial call");

    auto nA = gcmA->nA();
    auto nO = gcmA->gcmO->nA();
    auto nhc_ice = gcmA->nhc();

    // Convert TOPO arrays to 1-D zero-based indexing
    // ...on elevation grid
    blitz::TinyVector<int,2> const shape_E2(nhc_ice, nA);
    blitz::Array<double,2> fhcE2(reshape(topoA.fhc, shape_E2));
    blitz::Array<int,2> undericeE2(reshape(topoE.underice, shape_E2));
    blitz::Array<double,2>  elevE2(reshape(topoA.elevE, shape_E2));
    // ...on atmosphere grid
    auto foceanA(reshape1(topoA.focean));
    auto flakeA(reshape1(topoA.flake));
    auto fgrndA(reshape1(topoA.fgrnd));
    auto fgiceA(reshape1(topoA.fgice));
    auto zatmoA(reshape1(topoA.zatmo));

    // Read the original topo file [Ocean grid]
    NcIO ncio(topoO_fname, 'r');
    blitz::Array<double,2> foceanO2(nc_read_blitz(ncio.nc, "FOCEAN"));
    blitz::Array<double,2> flakeO2(nc_read_blitz(ncio.nc, "FLAKE"));
    blitz::Array<double,2> fgrndO2(nc_read_blitz(ncio.nc, "FGRND"));
    blitz::Array<double,2> fgiceO2(nc_read_blitz(ncio.nc, "FGICE"));
    blitz::Array<double,2> zatmoO2(nc_read_blitz(ncio.nc, "ZATMO"));

    blitz::Array<double,1> &foceanO(regrid1(foceanO2));
    blitz::Array<double,1> &flakeO(regrid1(flakeO2));
    blitz::Array<double,1> &fgrndO(regrid1(fgrndO2));
    blitz::Array<double,1> &fgiceO(regrid1(fgiceO2));
    blitz::Array<double,1> &zatmoO(regrid1(zatmoO2));
    ncio.close();

    // Keep track of which gridcells have been changed
    blitz::Array<char,1> changedO(nO);
    changedO = 0;

    // --------------------------------------
    // Add ice sheet to foceanO (and call it foceanOp)
    blitz::Array<double,1> &foceanOp(gcmA->foceanAOp);
    blitz::Array<double,1> &foceanOm(gcmA->foceanAOm);
    foceanOp = foceanO;
    update_foceanOp(gcmO, elevmasks, foceanOp, changedO);

    // --------------------------------------
    // Add ice to the surface type
    update_fgiceO(gcmc, gcmO, fgiceO, changedO);

    // --------------------------------------
    // Adjust fgrnd to make it all sum to 1; and round foceanOm at the same time
    for (int i=0; i<nO; ++i) {
        if (changedO(i)) {
            flakeO(i) = 0.;
            if (foceanOp(i) >= 0.5)
                foceanOm(i) = 1.;
                fgrndO(i) = 0.;
                fgiceO(i) = 0.;
            } else {
                foceanOm(i) = 0.;
                fgrnd(i) = 1. - fgiceO(i);
            }
        }
    }

    // Store the initial FOCEAN for ModelE, since it cannot change later.
    if (first) foceanOm0 = foceanOm;

    // ----------------------------------------------------------
    // ----------------------------------------------------------
    // Now we are ready to use regrid matrices

    // =====================================================
    // Regrid TOPO to Atmosphere grid
    HntrGrid const &hntrA(*cast_Grid_LonLat(&*this->gridA)->hntr);
    HntrGrid const &hntrO(*cast_Grid_LonLat(&*this->gcmO->gridA)->hntr);
    Hntr hntrAvO({&hntrA, &hntrO});
    TupleListT<double,2> AvO_tuples;
    hntrAvO.scaled_regrid(AvO_tuples);
    auto AvO_e(AvO_tuples.to_eigen());

    to_eigen(foceanA) = AvO_e * foceanOm;
    to_eigen(flakeA) = AvO_e * flakeO;
    to_eigen(fgrndA) = AvO_e * fgrnd;
    to_eigen(fgiceA) = AvO_e * fgice;
    to_eigen(zatmoA) = AvO_e * zatmo;

    // =====================================================
    // ---------- Compute elevE and AvE (aka fhc)

    // Whole-earth arrays, spanning all ice sheets
    TupleListT<double,2> AvE_global_tp;
    blitz::Array<double,1> elevE_global;

    // Compute elevE and AvE (aka fhc)
    for (auto sheet=sheets_index.begin(); sheet != sheets_index.end(); ++sheet) {
        // Construct an elevmaskI for ice sheet, =nan off ice sheet
        IceCoupler *icec(gcmc->ice_coupler(*sheet));
        blitz::Array elevmaskI(icec->nI());
        for (int iI=0; iI<icec->nI(); ++iI) {
            auto const m(icec->maskI(iI));
            if (m==MASK_GROUNDED_IE || m==MASK_FLOATING_ICE) {
                elevmaskI(i) = icec->elevI(iI);
            } else {
                elevmaskI(i) = nan;
            }
        }

        // Get regrid matrice needed to compute global stuff
        RegridMatrices rm(regrid_matrices(*sheet, elevmaskI));
        SparseSetT dimA, dimE, dimI;
        RegridMatrices::Params params;
            params.scale = true;
            params.correctA = false;
            params.sigma = ...;    // TODO: Set smoothing!
        auto AvI(rmA.matrix('AvI', {&dimA, &dimI}, params);
        auto EvI(rmA.matrix('EvI', {&dimE, &dimI}, params);
        auto AvE(rmA.matrix('AvE', {&dimE, &dimA}, params);

        // Merge local and global AvE
        spsparse::spcopy(
            spsparse::to_sparse({&dimA, &dimE},
            ref(AvE_global_tp)),
            AvE);


        // Merge local and global elevE
        auto elevE(EvI.apply(elevmaskI, nan, true));
        spsparse::spcopy(
            spsparse::to_sparse({&dimE},
            blitz_existing(elevE_global)),
            elevE);
    }

    // Create matrix that works directly on sparse-indexed vectors
    AvE_global_e = AvE_global_tp.to_eigen();
    EigenColVectorT elevA_global_e(AvE_global_e * eigen_col_vector(elevE_global));
    auto elevA_global(to_blitz(elevA_global_e));    // Sparse indexing

    // =======================================================
    // ----------- Set TOPO variables

    // =======================================================
    // ----------- Set up elevation class structure
    HCSegmentData const &legacy(get_segment(hc_segments, "legacy"));
    HCSegmentData const &sealand(get_segment(hc_segments, "sealand"));
    HCSegmentData const &ec(get_segment(hc_segments, "ec"));

    // Set up elevation class segments: fhc, underice, elevI
    fhcE2 = 0;
    undericeE2 = 0;
    elevE2 = 0;

    // ------- Segment 0: Legacy Segment
    for (int ihc=legacy.base; ihc<legacy.end; ++ihc) {
        for (int iA_s=0; iA_s<dimA.sparse_extent(); ++iA_s) {
            fhc(ihc,iA_s) = 1.;    // Legacy ice for Greenland and Antarctica.
            underice(ihc,iA_s) = UI_NOTHING;
        }

        for (int iA_s=0; iA_s<dimA.sparse_extent(); ++iA_s) {
            elevE(ihc,iA_s) = zatmoA(iA_s);
        }

        // overaly...
        for (int iA_d=0; iA_d<dimA.dense_extent(); ++iA_d) {
            int iA_s = dimA.to_sparse(iA_d);
            elevE(ihc,iA_s) = elevA(iA_s);
        }
    }}



    // ------- Segment 1: SeaLand Segment
    // ihc=0: Non-ice portion of grid cell at sea level

    // FHC is fraction of ICE-COVERED area in this elevation class
    // Therefore, FHC=0 for the sea portion of the SeaLand Segment
    // NOT: fhc[sealand.base,_maskA] = 1.-fgice[_maskA]
    // We could do away with this EC altogether because it's not used.
    for (int iA_d=0; iA_d<dimA.dense_extent(); ++iA_d) {
        int iA_s = dimA.to_sparse(iA_d);

        fhc(sealand.base, iA_s) = 0.
        underice(sealand.base, iA_s) = 0
        elevE(sealand.base, iA_s) = 0.
    };

    // ihc=1: Ice portion of grid cell at mean for the ice portion
    // FHC is fraction of ICE-COVERED area in this elevation class
    // Therefore, FHC=1 for the land portion of the SeaLand Segment
    // NOT: fhc[sealand.base+1,_maskA] = fgice[_maskA]
    for (int iA_d=0; iA_d<dimA.dense_extent(); ++iA_d) {
        int iA_s = dimA.to_sparse(iA_d);

        fhc(sealand.base+1, iA_s) = 1.
        underice(sealand.base+1, iA_s) = UI_NOTHING
        elevE(sealand.base+1, iA_s) = elevA[_maskA]
    }

    // ---------- Segment 2: Full Elevation Classes
    for (auto ii=AvE_global.begin(); ii != AvE_global.end(); ++ii) {
        auto const iA = ii->index(0);
        auto const iE = ii->index(1);

        auto const iE_tuple(indexingHC.index_to_tuple(iE));
        auto const iA2 = iE_tuple[0];
        auto const ihc = iE_tuple[1];

        if (iA2 != iA) (*icebin_error)(-1,
            "Matrix is non-local: iA=%d, iE=%d, iA2=%d", (int)iA, (int)iE, (int)iA2);

        fhc2(ec_base+ihc, iA) = ii->value();
        underice2(ec_base+ihc, iA) = UI_ICEBIN;
    }

    for (int ihc=0; ihc<elevE2.extent(0); ++ihc) {
    for (int iA=0; iA<elevE2.extent(1); ++iA) {
        elevE2(ec_base+ihc, iA) = hcdefs_ice(ihc);
    }}


    // ==================================================
    // 4) Fix bottom of atmosphere

    // Atmosphere sees mean of elevations over entire grid cell
    for (int iA_d=0; iA_d<dimA.dense_extent(); ++iA_d) {
        int iA_s = dimA.to_sparse(iA_d);

        zatmo_m(iA_s) = elevE2(legacy.base, iA_s);
    }







Input Files:
   1. TOPO output (O)             [Does not change]
   2. elevI / maskI (I)           [Set in GCMCoupler::couple()]
   3. IceBin INput File (O)

Initial Condition
-----------------
Inputs:
    Continent mask (cont_elevI)            [Use ice sheets elevI overlaid on top of original elevI0]
    TOPO output minus ice sheets (on O)

Outputs:
    TOPO outputs w/ ice sheets (on A)
        incl. underice
    GCMRegridder_ModelE (on A)
        --> fhc

Subsequent Coupling Steps
-------------------------
Inputs:
    Updated Continent mask (cont_elevI)
