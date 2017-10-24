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
blitz::Array<double,1> &foceanOp,    // 0-based array
std::map<std::string, blitz::Array<double,1>> const &cont_elevIs)
{
    auto nO(gcmO.nA());
    Grid const &gridO = *gcmO.gridA;

    // ------- Set up sO, scalaing factor for Ocean grid
    blitz::Array<double,1> sO(nO);    // 1. / |grid cells in O|
    sO = nan;
    for (auto cell=gridA->cells.begin(); cell != gridA->cells.end(); ++cell) {
        int iO_s = cell->index;
        sO(iO_s) = 1. / cell->native_area;
    }

    // --------------------- Compute fcontOp
    blitz::Array<double,1> fcontOp;

    for (auto sheet_name=sheets_index.begin(); sheet_name != sheets_index.end(); ++sheet_name) {
        // configure regridder for full continent (not just ice sheet)
        IceRegridder *iceO(gcmO.ice_regridder(*sheet_name));
        iceO->set_elevI(cont_elevIs.at(*sheet_name));

        // Get OvI for continental cells
        RegridMatrices rm(gcmO.regrid_matrices(*sheet_name));
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
        ibmisc::TmpAlloc tmp;
        blitz::Array<double,1> fcontOp_d(OvI.apply(fcontI_d, 0., true));    // force_conservation set to true by default, and it probably doesn't matter; but I think it should be false here.

        // Interpolate into foceanOp_s
        for (int iO_d=0; iO_d<fcontOp_d.extent(0); ++iO_d) {
            auto const iO_s = dimO.to_sparse(iO_d);
            foceanOp(iO_s) -= fcontOp_d(iO_d);
        }
    }
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
