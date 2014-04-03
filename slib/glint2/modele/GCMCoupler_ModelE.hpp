#pragma once

namespace glint2 {
namespace modele {

class GCMCoupler_ModelE : public GCMCoupler
{
	BOOST_ENUM_VALUES( CouplingType, int,
		/** GCM reports top T boundary condition to ice sheet.  This is
		always available. */
		(DIRICHLET_BC) (0)

		/** GCM reports energy fluxes at top of ice sheet.  This is only
		available on some ice models. */
		(NEUMANN_BC) (1)
	);

	GCMCoupler_ModelE(IceModel::GCMParams const &_gcm_params) :
		GCMCoupler(GCMCoupler::Type::MODELE, _gcm_params)
	{
		std::string descr;
		gcm_inputs.push_back(CoupledField("smb", "Surface mass balance", "kg m-2"));
		gcm_inputs.push_back(CoupledField("seb", "Latent heat flux, integrated", "J m-2"));
		gcm_inputs.push_back(CoupledField("tg2", "T of bottom layer of snow/firn, "C"));
		gcm_inputs_enum.set_names(coupled_field_names(gcm_inputs));
	}

	virtual void read_from_netcdf(
		NcFile &nc,
		std::string const &vname,
		std::vector<std::string> const &sheet_names,
	    giss::MapDict<std::string, IceSheet> &sheets);



};


}}

