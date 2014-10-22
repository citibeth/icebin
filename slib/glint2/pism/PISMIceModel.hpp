#pragma once

// --------------------------------
// PISM Includes... want to be included first
#include <petsc.h>
#include "IceGrid.hh"
#include "iceModel.hh"

#include "pism_options.hh"
#include "PAFactory.hh"
#include "POFactory.hh"
#include "PSFactory.hh"

#include <PISMTime.hh>
// --------------------------------
#include <boost/filesystem.hpp>
#include <glint2/pism/PSConstantGLINT2.hpp>
#include <glint2/pism/NullTransportHydrology.hpp>
#include <glint2/pism/MassEnergyBudget.hpp>

namespace glint2 {
namespace gpism {


/** This is the GLINT2 customized version of PISM's pism::IceModel class.

See https://github.com/pism/pism/issues/219

Here's a short term solution, though: create a new class
PISMIceModel derived from IceModel and re-implement
IceModel::allocate_couplers(). In it, set
IceModel::external_surface_model and IceModel::external_ocean_model as
you see fit (IceModel will not de-allocate a surface (ocean) model if
it is set to true) and allocate PSConstantGLINT2. You might also want
to add PISMIceModel::get_surface_model() which returns
IceModel::surface to get access to PSConstantGLINT2 from outside of
PISMIceModel.
*/

class PISMIceModel : public pism::IceModel
{
	friend class IceModel_PISM;
public:
	typedef pism::IceModel super;
	struct Params {
		double time_start_s;
		boost::filesystem::path output_dir;
	};
	Params const params;
protected:

	MassEnergyBudget base;		// Cumulative totals at start of this timestep
	MassEnergyBudget cur;		// Cumulative totals now
	MassEnergyBudget rate;		// At end of coupling timestep, set to (cur - base) / dt

	// Output variables prepared for return to GCM
	MassEnthVec2S basal_runoff;	// = melt_grounded + melt_floating [kg m-2 s-1]
	IceModelVec2S ice_surface_enth;		// Specific enthalpy of top surface of the ice [J kg-1]
	IceModelVec2S ice_surface_enth_depth;	// Depth below surface at which ice_surface_enth is recorded [m]
	


protected:
	// see iceModel.cc
	virtual PetscErrorCode createVecs();
	virtual PetscErrorCode allocate_internal_objects();

public:
	virtual PetscErrorCode massContExplicitStep();
	virtual PetscErrorCode accumulateFluxes_massContExplicitStep(
		int i, int j,
		double surface_mass_balance,		   // [m s-1] ice equivalent
		double meltrate_grounded,			  // [m s-1] ice equivalent
		double meltrate_floating,			  // [m s-1] ice equivalent
		double divQ_SIA,					   // [m s-1] ice equivalent
		double divQ_SSA,					   // [m s-1] ice equivalent
		double Href_to_H_flux,				 // [m s-1] ice equivalent
		double nonneg_rule_flux);			  // [m s-1] ice equivalent
private:
	// Temporary variables inside massContExplicitStep()
	double _ice_density;		// From config
	double _meter_per_s_to_kg_per_m2;


private:
	// Utility function
	PetscErrorCode prepare_nc(std::string const &fname, std::unique_ptr<pism::PIO> &nc);

public:

	/** @param t0 Time of last time we coupled. */
	PetscErrorCode set_rate(double dt);


	std::unique_ptr<pism::PIO> pre_mass_nc;	//!< Write variables every time massContPostHook() is called.
	std::unique_ptr<pism::PIO> post_mass_nc;
	std::unique_ptr<pism::PIO> pre_energy_nc;
	std::unique_ptr<pism::PIO> post_energy_nc;

	// see iceModel.cc for implementation of constructor and destructor:
	/** @param gcm_params Pointer to IceModel::gcm_params.  Lives at least as long as this object. */
	PISMIceModel(pism::IceGrid &g, pism::Config &config, pism::Config &overrides, Params const &params);
	virtual ~PISMIceModel(); // must be virtual merely because some members are virtual

	virtual PetscErrorCode allocate_enthalpy_converter();
	virtual PetscErrorCode allocate_subglacial_hydrology();
	virtual PetscErrorCode allocate_couplers();
	virtual PetscErrorCode grid_setup();
	virtual PetscErrorCode misc_setup();

	PetscErrorCode compute_enth2(pism::IceModelVec2S &enth2, pism::IceModelVec2S &mass2);

	/** @return Our instance of PSConstantGLINT2 */
	PSConstantGLINT2 *ps_constant_glint2()
		{ return dynamic_cast<PSConstantGLINT2 *>(surface); }
	NullTransportHydrology *null_hydrology()
		{ return dynamic_cast<NullTransportHydrology *>(pism::IceModel::subglacial_hydrology); }


	/** @return Current time for mass timestepping */
	double mass_t() const { return grid.time->current(); }
	/** @return Current time for enthalpy timestepping */
	double enthalpy_t() const { return t_TempAge; }

	PetscErrorCode massContPreHook();
	PetscErrorCode massContPostHook();
	// Pre and post for energy
	PetscErrorCode energyStep();

	PetscErrorCode write_post_energy(double time_s);

};

}}
