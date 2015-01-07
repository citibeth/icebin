#include <mpi.h>		// Must be first

#include <glint2/pism/GLINT2EnthalpyConverter.hpp>

#include <glint2/pism/IceModel_PISM.hpp>
#include <glint2/modele/GCMCoupler_ModelE.hpp>
#include <glint2/contracts/contracts.hpp>

using namespace giss;
using namespace glint2::modele;
using namespace glint2::gpism;

// --------------------------------------------------------
namespace glint2 {
namespace gpism {

static double const nan = std::numeric_limits<double>::quiet_NaN();

/** GCM-specific contract */
void IceModel_PISM::setup_contracts_modele()
{
	// Get arguments we need from coupler
	auto coupler(dynamic_cast<GCMCoupler_ModelE const *>(this->coupler));
	auto params(dynamic_cast<GCMPerIceSheetParams_ModelE const *>(this->gcm_per_ice_sheet_params.get()));


	printf("BEGIN IceModel_PISM::setup_contracts_modele\n");
	IceModel &model(*this);

	// =========== Transfer  constants
	transfer_constant("standard_gravity", "constant::grav");
	transfer_constant("beta_CC", "seaice::dtdp", -1.0);
	transfer_constant("water_melting_point_temperature", "constant::tf");
	transfer_constant("water_latent_heat_fusion", "constant::lhm");
	transfer_constant("water_specific_heat_capacity", "constant::shw");

	transfer_constant("ice_density", "constant::rhoi");
	transfer_constant("ice_thermal_conductivity", "seaice::alami0");
	transfer_constant("ice_specific_heat_capacity", "constant::shi");
	transfer_constant("fresh_water_density", "constant::rhow");
	transfer_constant("sea_water_density", "constant::rhows");
	transfer_constant("standard_gravity", "constant::grav");
	transfer_constant("ideal_gas_constant", "constant::gasc");

	Z2LI = coupler->gcm_constants.get_as("landice::z2li", "m");

	// To set this, see (in ModelE): Function SHCGS in ocnfuntab.f is
	// used for the Russell ocean. I. The simple models use SHW=4185.
	// This probably doesn't matter much at this point (May 2014)
	//	    transfer_constant("sea_water_specific_heat_capacity", "");

	/* The following constants were not transferred 
	pism_config:fill_value = -2e9;
	pism_config:fill_value_doc = "_FillValue used when saving diagnostic quantities";
	*/

	// In PISM and ModelE Clausius-Clapeyron equation, surfce_pressure is the DIFFERENCE
	// from 1atm.  Thus, surface_pressure=0 implies the ice sheet existing at 1atm
	set_constant("surface_pressure", 0, "Pa");		// Match ModelE thermodynam

	// No need to set enthalpy_reference_temperature.  pism::EnthalpyConverter is used (below)
	// to convert enthalpy values between ModelE and PISM.
	// transfer_constant("enthalpy_reference_temperature", "enthalpy_reference_temperature");


	// ============ GCM -> Ice
	CouplingContract &ice_input(contract[IceModel::INPUT]);

	// ------ Decide on the coupling contract for this ice sheet
	ice_input.add_field("smb_mass", nan, "kg m-2 s-1", contracts::ICE,
		"'Surface Mass Balance' over the coupling interval.\n"
		"Convention: Down is positive");
	ice_input.add_field("smb_enth", nan, "W m-2", contracts::ICE,
		"Advective enthalpy associated with smb_mass."
		"Convention: Down is positive");

	ice_input.add_field("surface_temp", nan, "K", contracts::ICE,
		"Mean temperature of the bottom of the ice surface model, over the "
		"coupling timestep.  This is for informational purposes ONLY, it is not "
		"used as a boundary condition.");

	ice_input.add_field("heat_flux", 0., "W m-2", contracts::ICE,
		"Conductive heat between ice sheet and snow/firn model on top of it.\n"
		"Convention: Down is positive");

	// Figure out the conversion between GCM and PISM enthalpy
	// ModelE's reference state is 1atm, 0C, 100% liquid water.
	// The enthalpy for that reference state would be the top end
	// of PISM's EnthalpyInterval.
	// NOTE: Pressure in PISM is RELATIVE to atmospheric pressure.
	//       Thus, p=0 is the correct to use at the top surface of
	//       the ice sheet (where ModelE operates).
	GLINT2EnthalpyConverter enth(*config);
	double const pressure = 0;
	double E_s, E_l;
	enth.getEnthalpyInterval(pressure, E_s, E_l);
	double const enth_modele_to_pism = E_l;		// (J/kg): Add to convert ModelE specific enthalpies (J/kg) to PISM specific enthalpies (J/kg)
	// NOTE: enth_modele_to_pism == 437000 J/kg
	if (pism_rank == 0) printf("enth_modele_to_pism = %g\n", enth_modele_to_pism);

	// ------------- Convert the contract to a var transformer
	VarTransformer &ice_input_vt(var_transformer[IceModel::INPUT]);
	ice_input_vt.set_names(VarTransformer::INPUTS, &coupler->gcm_outputs);
	ice_input_vt.set_names(VarTransformer::OUTPUTS, &ice_input);
	ice_input_vt.set_names(VarTransformer::SCALARS, &coupler->ice_input_scalars);
	ice_input_vt.allocate();

	// Add some recipes for gcm_to_ice
	std::string out;
	bool ok = true;
	ok = ok && ice_input_vt.set("smb_mass", "lismb", "unit", 1.0);

	// enthalpy flux (PISM) = liseb + enth_modele_to_pism * lismb
	ok = ok && ice_input_vt.set("smb_enth", "liseb", "unit", 1.0);
	ok = ok && ice_input_vt.set("smb_enth", "lismb", "unit", enth_modele_to_pism);

	ok = ok && ice_input_vt.set("surface_temp", "litg2", "unit", 1.0);
	ok = ok && ice_input_vt.set("surface_temp", "unit", "unit", C2K);	// +273.15

	ok = ok && ice_input_vt.set("heat_flux", "lif2", "unit", 1.0);

	// ============== Ice -> GCM
	CouplingContract &ice_output(contract[IceModel::OUTPUT]);

	// Glint2 requires that all ice models return elev2, so that it can regrid in the vertical.

	ice_output.add_field("usurf", "m", contracts::ICE|contracts::INITIAL, "ice upper surface elevation");	// See ice_surface_elevation in iceModel.cc

	ice_output.add_field("ice_surface_enth", "J kg-1", contracts::ICE|contracts::INITIAL, "");
	ice_output.add_field("ice_surface_enth_depth", "m", contracts::ICE|contracts::INITIAL, "");

	ice_output.add_field("effective_surface_temp", "degC", contracts::ICE, "");

	ice_output.add_field("basal_frictional_heating", "W m-2", contracts::ICE, "");
	ice_output.add_field("strain_heating", "W m-2", contracts::ICE, "");
	ice_output.add_field("geothermal_flux", "W m-2", contracts::ICE, "");
	ice_output.add_field("upward_geothermal_flux", "W m-2", contracts::ICE, "");

	ice_output.add_field("calving.mass", "kg m-2 s-1", contracts::ICE, "");
	ice_output.add_field("calving.enth", "W m-2", contracts::ICE, "");
	ice_output.add_field("glint2_smb.mass", "kg m-2 s-1", contracts::ICE, "");
	ice_output.add_field("glint2_smb.enth", "W m-2", contracts::ICE, "");
	ice_output.add_field("pism_smb.mass", "kg m-2 s-1", contracts::ICE, "");
	ice_output.add_field("pism_smb.enth", "W m-2", contracts::ICE, "");

	// basal_runoff (GCM input) = melt_grounded + melt_floatig (PISM outputs)
	ice_output.add_field("melt_grounded.mass", "kg m-2 s-1", contracts::ICE, "");
	ice_output.add_field("melt_grounded.enth", "W m-2", contracts::ICE, "");
	ice_output.add_field("melt_floating.mass", "kg m-2 s-1", contracts::ICE, "");
	ice_output.add_field("melt_floating.enth", "W m-2", contracts::ICE, "");

	ice_output.add_field("internal_advection.mass", "kg m-2 s-1", contracts::ICE, "");
	ice_output.add_field("internal_advection.enth", "W m-2", contracts::ICE, "");
	ice_output.add_field("epsilon.mass", "kg m-2 s-1", contracts::ICE, "");
	ice_output.add_field("epsilon.enth", "W m-2", contracts::ICE, "");

	ice_output.add_field("unit", "", 0, "Dimensionless identity");

	std::cout << "========= Ice Model Outputs (" << model.name << ") modele_pism.cpp:" << std::endl;
	std::cout << ice_output << std::endl;

	// ------- Variable and unit conversions, Ice -> GCM
	{VarTransformer &vt(var_transformer[IceModel::OUTPUT]);

	vt.set_names(VarTransformer::INPUTS, &ice_output);
	vt.set_names(VarTransformer::OUTPUTS, &coupler->gcm_inputs);
	vt.set_names(VarTransformer::SCALARS, &coupler->ice_input_scalars);
	vt.allocate();

//	ok = ok && vt.set("elev2", "usurf", "unit", 1.0);
	ok = ok && vt.set("elev1", "usurf", "unit", 1.0);

	// For specific enthalpy: Enth_e = Enth_p - enth_modele_to_pism
	// where X_e is ModelE and X_p is PISM
	ok = ok && vt.set("ice_surface_enth", "ice_surface_enth", "unit", 1.0);
	ok = ok && vt.set("ice_surface_enth", "unit", "unit", -enth_modele_to_pism);
	ok = ok && vt.set("ice_surface_enth_depth", "ice_surface_enth_depth", "unit", 1.0);



	ok = ok && vt.set("basal_frictional_heating", "basal_frictional_heating", "unit", 1.0);
	ok = ok && vt.set("strain_heating", "strain_heating", "unit", 1.0);

	ok = ok && vt.set("geothermal_flux", "geothermal_flux", "unit", 1.0);
	ok = ok && vt.set("upward_geothermal_flux", "upward_geothermal_flux", "unit", 1.0);

#if 0
	ok = ok && vt.set("glint2_smb.mass", "glint2_smb.mass", "unit", 1.0);
	ok = ok && vt.set("glint2_smb.enth", "glint2_smb.enth", "unit", 1.0);
	ok = ok && vt.set("glint2_smb.enth", "glint2_smb.mass", "unit", -enth_modele_to_pism);
#endif

	ok = ok && vt.set("basal_runoff.mass", "melt_grounded.mass", "unit", 1.0);
	ok = ok && vt.set("basal_runoff.enth", "melt_grounded.enth", "unit", 1.0);
	ok = ok && vt.set("basal_runoff.enth", "melt_grounded.mass", "unit", -enth_modele_to_pism);

	ok = ok && vt.set("basal_runoff.mass", "melt_floating.mass", "unit", 1.0);
	ok = ok && vt.set("basal_runoff.enth", "melt_floating.enth", "unit", 1.0);
	ok = ok && vt.set("basal_runoff.enth", "melt_floating.mass", "unit", -enth_modele_to_pism);

	ok = ok && vt.set("calving.mass", "calving.mass", "unit", 1.0);
	ok = ok && vt.set("calving.enth", "calving.enth", "unit", 1.0);
	ok = ok && vt.set("calving.enth", "calving.mass", "unit", -enth_modele_to_pism);


	ok = ok && vt.set("internal_advection.mass", "internal_advection.mass", "unit", 1.0);
	ok = ok && vt.set("internal_advection.enth", "internal_advection.enth", "unit", 1.0);
	ok = ok && vt.set("internal_advection.enth", "internal_advection.mass", "unit", -enth_modele_to_pism);


	ok = ok && vt.set("epsilon.mass", "epsilon.mass", "unit", 1.0);
	ok = ok && vt.set("epsilon.enth", "epsilon.enth", "unit", 1.0);
	ok = ok && vt.set("epsilon.enth", "epsilon.mass", "unit", -enth_modele_to_pism);

	}

	// Catch all our errors at once
	if (!ok) throw std::exception();
	printf("END IceModel_PISM::setup_contracts_modele\n");
}


}}		// namespace glint2::gpism
// --------------------------------------------------------


// Contracts should also specify how constants are agreed upon between the two parties.
// 
// PISM needs at least the following constants:
// 
// EnthalpyConverter::EnthalpyConverter(const PISMConfig &config) {
//   beta  = config.get("beta_CC");                                 // K Pa-1
//   c_i   = config.get("ice_specific_heat_capacity");              // J kg-1 K-1
//   g     = config.get("standard_gravity");                        // m s-2
//   L     = config.get("water_latent_heat_fusion");                // J kg-1
//   p_air = config.get("surface_pressure");                        // Pa
//   rho_i = config.get("ice_density");                             // kg m-3
//   T_melting = config.get("water_melting_point_temperature");       // K  
//   T_tol = config.get("cold_mode_is_temperate_ice_tolerance");    // K 
//   T_0   = config.get("enthalpy_converter_reference_temperature");// K  
// 
//   do_cold_ice_methods  = config.get_flag("do_cold_ice_methods");
// }
