#include <glint2/modele/GCMCoupler_ModelE.hpp>

namespace glint2 {
namespace modele {

void GCMCoupler_ModelE::read_from_netcdf(
	NcFile &nc,
	std::string const &vname,
	std::vector<std::string> const &sheet_names,
    giss::MapDict<std::string, IceSheet> &sheets)
{
	printf("BEGIN GCMCoupler_ModelE::read_from_netcdf()\n");

	// Fills in models
	GCMCoupler::read_from_netcdf();

	// Read GCM-specific coupling parameters
	// Set the contract for each ice sheet, based on:
	//   (a) GCM-specific coupling parameters (to be read),
	//   (b) The type of ice model
	int i = 0;
	for (auto name = sheet_names.begin(); name != sheet_names.end(); ++name, ++i) {
		PerSheet &sheet = per_sheet[i];

		std::string sheet_vname = vname + "." + *name;
		auto gcm_var = nc.get_var((sheet_vname + ".modele").c_str());

		CouplingType coupling_type = giss::parse_enum<CouplingType>(
			giss::get_att(gcm_var, "coupling_type")->as_string(0));

		// --------------------------------------------------------------
		// ------ Decide on the coupling contract for this ice sheet
		CouplingContract &con(sheet.contract);

		// GCM to Ice
		cf = giss::get_cfname("land_ice_surface_specific_mass_balance_flux", "kg m-2 s-1");
		con.gcm_to_ice.push_back(CoupledField(cf));
		cf = giss::get_cfname("surface_downward_latent_heat_flux", "W m-2");
		con.gcm_to_ice.push_back(CoupledField(cf));
		switch(coupling_type.index()) {
			case CouplingType::DIRICHLET :
				cf = giss::get_cfname("surface_temperature", "K");
			break;
			case CouplingType::NEUMANN :
				cf = giss::get_cfname("surface_downward_sensible_heat_flux", "W m-2");
			break;
		}
		con.gcm_to_ice.push_back(CoupledField(cf));

		// Ice to GCM
		std::string descr;
		con.ice_to_gcm.push_back(CoupledField("upward_geothermal_flux_sum", descr,	"J m-2"));
		con.ice_to_gcm.push_back(CoupledField("geothermal_flux_sum", descr,			"J m-2"));
		con.ice_to_gcm.push_back(CoupledField("basal_frictional_heating_sum", descr,"J m-2"));
		con.ice_to_gcm.push_back(CoupledField("strain_heating_sum",descr,			"J m-2"));
		con.ice_to_gcm.push_back(CoupledField("total_enthalpy",descr,				"J m-2"));
		// --------------------------------------------------------------

		// ------------- Convert the contract to a var transformer
		VarTransformer &vt(sheet.vt_gcm_to_ice);
		std::vector<string> names;
		vt.set_names(VarTransformer::OUTPUTS, coupled_field_names(con.gcm_to_ice));
		names = coupled_field_names(gcm_inputs);
		names.push_back("unit");
		vt.set_names(VarTransformer::INPUTS, std::move(names));
		vt.set_names(VarTransformer::SCALARS, {
			"by_dt",		// 1.0 / (Coupling timestep)
			"unit"});

		// Add some recipes for gcm_to_ice
		std::string out;
		out = "land_ice_surface_specific_mass_balance_flux";
			vt.set(out, "smb", "by_dt", 1.0);
		out = "surface_downward_latent_heat_flux";
			vt.set(out, "seb", "by_dt", 1.0);
		out = "surface_temperature";	// K
			vt.set(out, "tg2", "unit", 1.0);
			vt.set(out, "unit", "unit", C2K);	// +273.15
		out = "surface_downward_sensible_heat_flux";	// W m-2
			// Zero for now

		// Print it out, see if the formulas make sense!
		std::cout << vt;

	}

	// Print out the contract
	i = 0;
	for (auto name = sheet_names.begin(); name != sheet_names.end(); ++name, ++i) {
		CouplingContract *con = contracts[i];
		std::out << "========= Contract for " << *name << std::endl;
		std::out << *con;
	}

	printf("END GCMCoupler_ModelE::read_from_netcdf()\n");
}

}}
