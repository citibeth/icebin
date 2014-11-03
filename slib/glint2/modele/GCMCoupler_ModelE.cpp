#include <mpi.h>		// Intel MPI wants to be first
#include <giss/memory.hpp>
#include <glint2/modele/GCMCoupler_ModelE.hpp>

namespace glint2 {
namespace modele {


GCMCoupler_ModelE::GCMCoupler_ModelE() :
	GCMCoupler(GCMCoupler::Type::MODELE)
{

	// ------------ GCM Outputs
	// The GCM must produce the same set of outputs, no matter what
	// ice model is being used
	gcm_outputs.add_field("lismb", "kg m-2 s-1", "Surface mass balance");
	gcm_outputs.add_field("liseb", "J m-2 s-1", "Ethalpy of SMB");
	// Specific enthalpy of the SMB computed as liseb / lismb (J kg-1)

	gcm_outputs.add_field("litg2", "degC", "T of bottom layer of snow/firn");
	gcm_outputs.add_field("unit", "", "Dimensionless identity");


	// ------------------------- GCM Inputs
	// The GCM must accept the same set of inputs, no matter what
	// ice moel is being used

	gcm_inputs.add_field("elev1", "m", "ATMOSPHERE", "ice upper surface elevation");

	// Enthalpy (T) of top of dynamic ice sheet used to compute heat
	// fluxes between ice surface model and the dynamic ice model (eg PISM)
	gcm_inputs.add_field("ice_surface_enth", "J kg-1", "ELEVATION", "");
	gcm_inputs.add_field("ice_surface_enth_depth", "m", "ELEVATION", "");

	// Mass/Energy outputs the GCM must deal with
	gcm_inputs.add_field("basal_runoff.mass", "kg m-2 s-1", "ATMOSPHERE", "");		// melt_grounded + melt_floating
	gcm_inputs.add_field("basal_runoff.enth", "W m-2", "ATMOSPHERE", "")
	gcm_inputs.add_field("calving.mass", "kg m-2 s-1", "ATMOSPHERE", "");
	gcm_inputs.add_field("calving.enth", "W m-2", "ATMOSPHERE", "");
	gcm_inputs.add_field("strain_heating", "W m-2", "ATMOSPHERE", "");

	// Unexplained mass/energy outputs, used to close the energy budget.
	gcm_inputs.add_field("epsilon.mass", "kg m-2 s-1", "ATMOSPHERE", "");
	gcm_inputs.add_field("epsilon.enth", "W m-2", "ATMOSPHERE", "");


	// ----------------- Scalars provided by the GCM
	// Scalars are things that can only be computed at the last minute
	// (eg, dt for a particular coupling timestep).  Constants that
	// can be computed at or before contract initialization time can
	// be placed directly into the VarTransformer.

// We don't need this, GCM is converting to s-1 on its own.
//	ice_input_scalars.add_field("by_dt", "s-1", "Inverse of coupling timestep");

	ice_input_scalars.add_field("unit", "", "Dimensionless identity");
//	gcm_input_scalars.add_field("unit", "", "Dimensionless identity");


#if 0
	ice_output_scalars.add_field("dt", "s", "Inverse of coupling timestep");
	ice_output_scalars.add_field("unit", "", "Dimensionless identity");
#endif

}


std::unique_ptr<GCMPerIceSheetParams>
GCMCoupler_ModelE::read_gcm_per_ice_sheet_params(
	NcFile &nc,
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
