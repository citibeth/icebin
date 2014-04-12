#pragma once

#include <glint2/GCMCoupler.hpp>

namespace glint2 {
namespace modele {


BOOST_ENUM_VALUES( ModelE_CouplingType, int,
	/** GCM reports top T boundary condition to ice sheet.  This is
	always available. */
	(DIRICHLET_BC) (0)

	/** GCM reports energy fluxes at top of ice sheet.  This is only
	available on some ice models. */
	(NEUMANN_BC) (1)
);

struct ContractParams_ModelE {
	ModelE_CouplingType coupling_type;
};

class GCMCoupler_ModelE : public GCMCoupler
{
public:
	/** Names of items used in the SCALARS dimension of VarTranslator */
	CouplingContract ice_input_scalars;

	GCMCoupler_ModelE(IceModel::GCMParams const &_gcm_params);

	virtual void read_from_netcdf(
		NcFile &nc,
		std::string const &vname,
		std::vector<std::string> const &sheet_names,
	    giss::MapDict<std::string, IceSheet> &sheets);

	virtual void setup_contracts(
		IceModel &mod,
		NcFile &nc,
		std::string const &sheet_vname);


};


}}

