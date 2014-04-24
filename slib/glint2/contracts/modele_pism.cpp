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
	ice_input.add_cfname("land_ice_surface_specific_mass_balance_flux", "kg m-2 s-1");
	ice_input.add_field("land_ice_surface_specific_mass_balance_water_fraction", "1", "");
	switch(params.coupling_type.index()) {
		case ModelE_CouplingType::DIRICHLET_BC :
			ice_input.add_cfname("surface_temperature", "K");
		break;
		case ModelE_CouplingType::NEUMANN_BC :
			ice_input.add_cfname("land_ice_surface_downward_conductive_heat_flux", "W m-2");
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
