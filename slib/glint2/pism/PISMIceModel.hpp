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
#include <glint2/pism/PSConstantGLINT2.hpp>
#include <glint2/pism/NullTransportHydrology.hpp>

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
	};
	Params const params;
protected:

	// Temporary stuff to hold returns from PISM
//	pism::IceModelVec2S strain_heating2;		//!< Rate of strain heating (W/m^2)
	pism::IceModelVec2S upward_geothermal_flux;	//!< Per-timestep upward geothermal flux (W/m^2)

	pism::IceModelVec2S basal_frictional_heating_sum;	//!< Total amount of basal friction heating (J/m^2)
	pism::IceModelVec2S strain_heating_sum;	//!< Total amount of strain heating (J/m^2)
	pism::IceModelVec2S geothermal_flux_sum;	//!< Total amount of geothermal energy (J/m^2)
	pism::IceModelVec2S upward_geothermal_flux_sum;	//!< Total amount of geothermal energy (J/m^2)
	pism::IceModelVec2S total_enthalpy;		//!< Total enthalpy of ice sheet (J/m^2)

protected:
	// see iceModel.cc
	virtual PetscErrorCode createVecs();

private:
	// Utility function
	PetscErrorCode prepare_nc(std::string const &fname, std::unique_ptr<pism::PIO> &nc);

public:

	std::unique_ptr<pism::PIO> pre_mass_nc;	//!< Write variables every time massContPostHook() is called.
	std::unique_ptr<pism::PIO> post_mass_nc;
	std::unique_ptr<pism::PIO> pre_energy_nc;
	std::unique_ptr<pism::PIO> post_energy_nc;

	// see iceModel.cc for implementation of constructor and destructor:
	/** @param gcm_params Pointer to IceModel::gcm_params.  Lives at least as long as this object. */
	PISMIceModel(pism::IceGrid &g, pism::PISMConfig &config, pism::PISMConfig &overrides, Params const &params);
	virtual ~PISMIceModel(); // must be virtual merely because some members are virtual

	virtual PetscErrorCode allocate_subglacial_hydrology();
	virtual PetscErrorCode allocate_couplers();
	virtual PetscErrorCode grid_setup();
    virtual PetscErrorCode misc_setup();

	virtual PetscErrorCode run_to(double time);

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
