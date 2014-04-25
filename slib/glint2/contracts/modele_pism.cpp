#include <mpi.h>		// Must be first

#include <glint2/pism/IceModel_PISM.hpp>
#include <glint2/modele/GCMCoupler_ModelE.hpp>

using namespace giss;
using namespace glint2::modele;
using namespace glint2::gpism;

// --------------------------------------------------------
namespace glint2 {
namespace gpism {

/** GCM-specific contract */
void IceModel_PISM::setup_contract_modele(
	glint2::modele::GCMCoupler_ModelE const &coupler,
	glint2::modele::ContractParams_ModelE const &params)
{
	printf("BEGIN IceModel_PISM::setup_contract_modele\n");
	IceModel &model(*this);

	// ============ GCM -> Ice
	CouplingContract &ice_input(contract[IceModel::INPUT]);

	// ------ Decide on the coupling contract for this ice sheet
	ice_input.add_field("land_ice_surface_downward_mass_flux", "kg m-2 s-1",
		"'Surface Mass Balance' over the coupling interval.\n"
		"Convention: Down is positive");
	ice_input.add_field("land_ice_surface_downward_enthalpy_flux", "J m-2",
		"Advective heat associated with land_ice_surface_downward_mass_flux."
		"Convention: Down is positive");
* Let's figure out for sure how to convert between ModelE and PISM Enthalpy
	switch(params.coupling_type.index()) {
		case ModelE_CouplingType::DIRICHLET_BC :
			ice_input.add_cfname("surface_temperature", "K");
		break;
		case ModelE_CouplingType::NEUMANN_BC :
			ice_input.add_field("land_ice_surface_downward_conductive_heat_flux", "W m-2",
				"Conductive heat between ice sheet and snow/firn model on top of it.\n"
				"Convention: Down is positive");
		break;
	}

	// ------------- Convert the contract to a var transformer
	VarTransformer &ice_input_vt(var_transformer[IceModel::INPUT]);
	ice_input_vt.set_names(VarTransformer::INPUTS, &coupler.gcm_outputs);
	ice_input_vt.set_names(VarTransformer::OUTPUTS, &ice_input);
	ice_input_vt.set_names(VarTransformer::SCALARS, &coupler.ice_input_scalars);
	ice_input_vt.allocate();

	// Add some recipes for gcm_to_ice
	std::string out;
	out = "land_ice_surface_specific_mass_balance_flux";
		ice_input_vt.set(out, "lismb", "by_dt", 1.0);
	out = "land_ice_surface_downward_advective_heat_flux";
		ice_input_vt.set(out, "liseb", "by_dt", 1.0);
	out = "surface_temperature";	// K
		ice_input_vt.set(out, "litg2", "by_dt", 1.0);
		ice_input_vt.set(out, "unit", "unit", C2K);	// +273.15
	out = "land_ice_surface_downward_conductive_heat_flux";	// W m-2
		// Zero for now

	// ============== Ice -> GCM
	CouplingContract &ice_output(contract[IceModel::OUTPUT]);
	ice_output.add_field("upward_geothermal_flux_sum", "J m-2", "");
	ice_output.add_field("geothermal_flux_sum", "J m-2", "");
	ice_output.add_field("basal_frictional_heating_sum", "J m-2", "");
	ice_output.add_field("strain_heating_sum", "J m-2", "");
	ice_output.add_field("total_enthalpy", "J m-2", "");
	ice_output.add_field("unit", "", "");

	// Outputs (Ice -> GCM) are same fields as inputs
	CouplingContract *gcm_inputs = new_CouplingContract();
	for (auto ii = ice_output.begin(); ii != ice_output.end(); ++ii) {
		gcm_inputs->add_field(*ii);
	}

	CouplingContract *ice_output_scalars = new_CouplingContract();
	ice_output_scalars->add_field("unit", "", "");

	VarTransformer &ice_output_vt(var_transformer[OUTPUT]);
	ice_output_vt.set_names(VarTransformer::INPUTS, &ice_output);
	ice_output_vt.set_names(VarTransformer::OUTPUTS, gcm_inputs);
	ice_output_vt.set_names(VarTransformer::SCALARS, ice_output_scalars);
	ice_output_vt.allocate();

	// Set up transformations: just copy inputs to outputs
	for (auto ii = ice_output.begin(); ii != ice_output.end(); ++ii) {
		ice_output_vt.set(ii->name, ii->name, "unit", 1.0);
	}

	// Now give our contracts to our dismal slave IceModel
	if (dismal.get()) {
		dismal->contract = contract;
	}

	printf("END IceModel_PISM::setup_contract_modele\n");
}

}}		// namespace glint2::gpism
// --------------------------------------------------------

Contracts should also specify how constants are agreed upon between the two parties.

PISM needs at least the following constants:

EnthalpyConverter::EnthalpyConverter(const PISMConfig &config) {
  beta  = config.get("beta_CC");                                 // K Pa-1
  c_i   = config.get("ice_specific_heat_capacity");              // J kg-1 K-1
  g     = config.get("standard_gravity");                        // m s-2
  L     = config.get("water_latent_heat_fusion");                // J kg-1
  p_air = config.get("surface_pressure");                        // Pa
  rho_i = config.get("ice_density");                             // kg m-3
  T_melting = config.get("water_melting_point_temperature");       // K  
  T_tol = config.get("cold_mode_is_temperate_ice_tolerance");    // K 
  T_0   = config.get("enthalpy_converter_reference_temperature");// K  

  do_cold_ice_methods  = config.get_flag("do_cold_ice_methods");
}