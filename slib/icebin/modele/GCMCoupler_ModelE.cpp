#include <mpi.h>        // Intel MPI wants to be first
#include <giss/memory.hpp>
#include <icebin/modele/GCMCoupler_ModelE.hpp>
#include <icebin/contracts/contracts.hpp>

namespace icebin {
namespace modele {

static double const nan = std::numeric_limits<double>::quiet_NaN();
using namespace icebin::contracts;

GCMCoupler_ModelE::GCMCoupler_ModelE() :
    GCMCoupler(GCMCoupler::Type::MODELE)
{

    // ------------ GCM Outputs
    // The GCM must produce the same set of outputs, no matter what
    // ice model is being used
    gcm_outputs.add_field("runo", nan, "kg m-2", ELEVATION,
        "Downward water flux through bottom layer");
    gcm_outputs.add_field("eruno", nan, "J m-2", ELEVATION,
        "Enthalpy of downward water flux through bottom layer");
#ifdef TRACERS_WATER
    gcm_outputs.add_field("trruno")
#endif
    gcm_outputs.add_field("deltah", nan, "J m-2", ELEVATION,
        "Enthalpy change of 'borrowed' layer");
    gcm_outputs.add_field("massxfer", nan, "kg m-2", ELEVATION,
        "Mass of ice being transferred Stieglitz --> Icebin");
    gcm_outputs.add_field("enthxfer", nan, "J m-2", ELEVATION,
        "Enthlpy of ice being transferred Stieglitz --> Icebin");
#ifdef TRACERS_WATER
    gcm_outputs.add_field("trxfer");
#endif
    gcm_outputs.add_field("volxfer", nan, "m^3 m-2", ELEVATION,
        "Volume of ice being transferred Stieglitz --> Icebin");

    gcm_outputs.add_field("unit", nan, "", 0, "Dimensionless identity");


    // ------------------------- GCM Inputs
    // ModelE sets this, via repeated calls to add_gcm_input_ij()
    // and add_gcm_input_ijhc().  See alloc_landic_com() in LANDICE_COM.f

    // ----------------- Scalars provided by the GCM
    // Scalars are things that can only be computed at the last minute
    // (eg, dt for a particular coupling timestep).  Constants that
    // can be computed at or before contract initialization time can
    // be placed directly into the VarTransformer.

    ice_input_scalars.add_field("by_dt", nan, "s-1", 1., "Inverse of coupling timestep");

    ice_input_scalars.add_field("unit", nan, "", 0, "Dimensionless identity");
//  gcm_input_scalars.add_field("unit", nan, "", 0, "Dimensionless identity");

}


std::unique_ptr<GCMPerIceSheetParams>
GCMCoupler_ModelE::read_gcm_per_ice_sheet_params(
    ibmisc::NcIO &ncio,
    std::string const &sheet_vname)
{


    // Read GCM-specific coupling parameters
    // Set the contract for each ice sheet, based on:
    //   (a) GCM-specific coupling parameters (to be read),
    //   (b) The type of ice model

    auto gcm_var = giss::get_var_safe(nc, (sheet_vname + ".modele").c_str());

    std::unique_ptr<GCMPerIceSheetParams_ModelE> params(
        new GCMPerIceSheetParams_ModelE());

    params->coupling_type = giss::parse_enum<ModelE_CouplingType>(
        giss::get_att(gcm_var, "coupling_type")->as_string(0));

    return giss::static_cast_unique_ptr<GCMPerIceSheetParams>(params);
}



void GCMCoupler_ModelE::setup_contracts(IceModel &ice_model) const
    { ice_model.setup_contracts_modele(); }


}}
