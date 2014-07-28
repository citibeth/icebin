#include <mpi.h>		// Intel MPI wants to be first
#include <giss/memory.hpp>
#include <glint2/modele/GCMCoupler_ModelE.hpp>

namespace glint2 {
namespace modele {


GCMCoupler_ModelE::GCMCoupler_ModelE() :
	GCMCoupler(GCMCoupler::Type::MODELE)
{
	gcm_outputs.add_field("lismb", "kg m-2 s-1", "Surface mass balance");
//	gcm_outputs.add_field("liseb", "J kg-1", "Specific ethalpy of SMB");
	gcm_outputs.add_field("liseb", "J m-2 s-1", "Ethalpy of SMB");
	gcm_outputs.add_field("litg2", "degC", "T of bottom layer of snow/firn");
	gcm_outputs.add_field("unit", "", "Dimensionless identity");

	ice_input_scalars.add_field("by_dt", "s-1", "Inverse of coupling timestep");
	ice_input_scalars.add_field("unit", "", "Dimensionless identity");

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
