#include <mpi.h>		// Must be first

#include <glint2/IceModel_DISMAL.hpp>
#include <glint2/modele/GCMCoupler_ModelE.hpp>

using namespace giss;
using namespace glint2::modele;

// --------------------------------------------------------
namespace glint2 {

/** GCM-specific contract */
void IceModel_DISMAL::setup_contracts_modele()
{
	// Get arguments we need from coupler
	auto coupler(dynamic_cast<GCMCoupler_ModelE const *>(this->coupler));
	auto params(dynamic_cast<GCMPerIceSheetParams_ModelE const *>(this->gcm_per_ice_sheet_params.get()));

	printf("BEGIN IceModel_DISMAL::setup_contract_modele\n");
	IceModel &model(*this);

	// Don't bother transferring any constants...

	// ============ GCM -> Ice
	CouplingContract &ice_input(contract[IceModel::INPUT]);

	std::string const MASS_FLUX = "surface_downward_mass_flux";
	std::string const ENTHALPY_FLUX = "surface_downward_enthalpy_flux";
	std::string const T = "surface_temperature";
	std::string const HEAT_FLUX = "surface_downward_conductive_heat_flux";

	// ------ Decide on the coupling contract for this ice sheet
	ice_input.add_field(MASS_FLUX, "kg m-2 s-1");
	ice_input.add_field(ENTHALPY_FLUX, "W m-2");
	switch(params->coupling_type.index()) {
		case ModelE_CouplingType::DIRICHLET_BC :
			ice_input.add_field(T, "K");
		break;
		case ModelE_CouplingType::NEUMANN_BC :
			ice_input.add_field(HEAT_FLUX, "W m-2");
		break;
	}

	// ------------- Convert the contract to a var transformer
	VarTransformer &ice_input_vt(var_transformer[IceModel::INPUT]);
	ice_input_vt.set_names(VarTransformer::INPUTS, &coupler->gcm_outputs);
	ice_input_vt.set_names(VarTransformer::OUTPUTS, &ice_input);
	ice_input_vt.set_names(VarTransformer::SCALARS, &coupler->ice_input_scalars);
	ice_input_vt.allocate();

	// Add some recipes for gcm_to_ice
	double const enth_modele_to_pism = 437000;		// Taken from run with PISM
	std::string out;
	bool ok = true;
	ok = ok && ice_input_vt.set(MASS_FLUX, "lismb", "unit", 1.0);

	ok = ok && ice_input_vt.set(ENTHALPY_FLUX, "liseb", "unit", 1.0);
	ok = ok && ice_input_vt.set(ENTHALPY_FLUX, "lismb", "unit", enth_modele_to_pism);

	switch(params->coupling_type.index()) {
		case ModelE_CouplingType::DIRICHLET_BC :
			ok = ok && ice_input_vt.set(T, "litg2", "unit", 1.0);
			ok = ok && ice_input_vt.set(T, "unit", "unit", C2K);	// +273.15
		break;
		case ModelE_CouplingType::NEUMANN_BC :
// Nothing for now...
//			ok = ok && ice_input_vt.set(HEAT_FLUX, "liseb", "unit", 1.0);
		break;
	}


	// ============== Ice -> GCM (none)
	CouplingContract &ice_output(contract[IceModel::OUTPUT]);

	// Outputs (Ice -> GCM) are same fields as inputs
	CouplingContract *gcm_inputs = new_CouplingContract();
	for (auto ii = ice_output.begin(); ii != ice_output.end(); ++ii) {
		gcm_inputs->add_field(*ii);
	}

	CouplingContract *ice_output_scalars = new_CouplingContract();
	ice_output_scalars->add_field("unit", "", 0, "");

	VarTransformer &ice_output_vt(var_transformer[OUTPUT]);
	ice_output_vt.set_names(VarTransformer::INPUTS, &ice_output);
	ice_output_vt.set_names(VarTransformer::OUTPUTS, gcm_inputs);
	ice_output_vt.set_names(VarTransformer::SCALARS, ice_output_scalars);
	ice_output_vt.allocate();

	// Set up transformations: just copy inputs to outputs
	for (auto ii = ice_output.begin(); ii != ice_output.end(); ++ii) {
		ok = ok && ice_output_vt.set(ii->name, ii->name, "unit", 1.0);
	}

	if (!ok) {
		printf("modele_dismal.cpp quitting due to errors.\n");
		throw std::exception();
	}

	printf("END IceModel_DISMAL::setup_contract_modele\n");
}

}		// namespace glint2
// --------------------------------------------------------
