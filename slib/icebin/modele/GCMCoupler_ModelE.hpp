#pragma once

#include <icebin/GCMCoupler.hpp>

namespace icebin {
namespace modele {


BOOST_ENUM_VALUES( ModelE_CouplingType, int,
	/** GCM reports top T boundary condition to ice sheet.  This is
	always available. */
	(DIRICHLET_BC) (0)

	/** GCM reports energy fluxes at top of ice sheet.  This is only
	available on some ice models. */
	(NEUMANN_BC) (1)
);


class GCMPerIceSheetParams_ModelE : public icebin::GCMPerIceSheetParams {
public:
	ModelE_CouplingType coupling_type;
};

class GCMCoupler_ModelE : public GCMCoupler
{
public:
	GCMCoupler_ModelE();

	/** Read per-ice-sheet parameters that depend on the type of GCMCoupler. */
	std::unique_ptr<GCMPerIceSheetParams>
	read_gcm_per_ice_sheet_params(
	    ibmisc::NcIO &ncio,
		std::string const &sheet_vname);

	/** Does contract setup for ONE IceModel instance.
	Calls throught to IceModel::setup_contract_xxx() */
	virtual void setup_contracts(IceModel &mod) const;

};


}}

