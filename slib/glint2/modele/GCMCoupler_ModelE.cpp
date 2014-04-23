#include <mpi.h>		// Intel MPI wants to be first
#include <glint2/modele/GCMCoupler_ModelE.hpp>

namespace glint2 {
namespace modele {


GCMCoupler_ModelE::GCMCoupler_ModelE(GCMParams const &_gcm_params) :
	GCMCoupler(GCMCoupler::Type::MODELE, _gcm_params)
{
	gcm_outputs.add_field("lismb", "kg m-2 s-1", "Surface mass balance");
	gcm_outputs.add_field("liseb", "W m-2", "Latent heat flux");
	gcm_outputs.add_field("litg2", "degC", "T of bottom layer of snow/firn");
	gcm_outputs.add_field("unit", "", "Dimensionless identity");

	ice_input_scalars.add_field("by_dt", "s-1", "Inverse of coupling timestep");
	ice_input_scalars.add_field("unit", "", "Dimensionless identity");

#if 0
	ice_output_scalars.add_field("dt", "s", "Inverse of coupling timestep");
	ice_output_scalars.add_field("unit", "", "Dimensionless identity");
#endif

}





void GCMCoupler_ModelE::setup_contracts(
	IceModel &model,
	NcFile &nc,
	std::string const &sheet_vname)
{
	printf("BEGIN GCMCoupler_ModelE::setup_contracts()\n");


	// Read GCM-specific coupling parameters
	// Set the contract for each ice sheet, based on:
	//   (a) GCM-specific coupling parameters (to be read),
	//   (b) The type of ice model

	auto gcm_var = giss::get_var_safe(nc, (sheet_vname + ".modele").c_str());

	ContractParams_ModelE params;
	params.coupling_type = giss::parse_enum<ModelE_CouplingType>(
		giss::get_att(gcm_var, "coupling_type")->as_string(0));

	model.setup_contract_modele(*this, params);
	model.finish_contract_setup();

	printf("END GCMCoupler_ModelE::setup_contracts()\n");
}


}}
