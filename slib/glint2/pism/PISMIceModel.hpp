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
#include <glint2/pism/VecBundleWriter.hpp>

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
	// (relevant ice model state to be exported)
	pism::IceModelVec2S M1,M2;
	pism::IceModelVec2S H1,H2;
	pism::IceModelVec2S V1,V2;

protected:
	// see iceModel.cc
	virtual PetscErrorCode createVecs();
//	virtual PetscErrorCode allocate_internal_objects();

public:
	virtual PetscErrorCode massContExplicitStep();
	virtual PetscErrorCode accumulateFluxes_massContExplicitStep(
		int i, int j,
		double surface_mass_balance,		   // [m s-1] ice equivalent (from PISM)
		double meltrate_grounded,			  // [m s-1] ice equivalent
		double meltrate_floating,			  // [m s-1] ice equivalent
		double divQ_SIA,					   // [m s-1] ice equivalent
		double divQ_SSA,					   // [m s-1] ice equivalent
		double Href_to_H_flux,				 // [m s-1] ice equivalent
		double nonneg_rule_flux);			  // [m s-1] ice equivalent
private:
	// Temporary variables inside massContExplicitStep()
	double _ice_density;		// From config
	double _meter_per_s_to_kg_per_m2;		// Conversion factor computed from _ice_density


private:
	// Utility function
//	PetscErrorCode prepare_nc(std::string const &fname, std::unique_ptr<pism::PIO> &nc);

public:

	/** @param t0 Time of last time we coupled. */
	PetscErrorCode set_rate(double dt);

	PetscErrorCode reset_rate();

	/** Various internal PISM variables at the point that PISM output
	is sent to the GCM */
	std::unique_ptr<VecBundleWriter> pism_out_state_nc;

	/** Various internal PISM variables just after input from the GCM
	is absorbed into state (but before running the PISM timestep) */
	std::unique_ptr<VecBundleWriter> pism_in_state_nc;


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

	PetscErrorCode prepare_outputs(double time_s);

	/** Read things out of the ice model that will be sent back BEFORE
	the first coupling timestep (eg, ice surface enthalpy) */
	PetscErrorCode prepare_initial_outputs();

	/** Merges surface temperature derived from Enth3 into any NaN values
	in the vector provided.
	@param deltah IN: Input from Glint2 (change in enthalpy of each grid
		cell over the timestep) [W m-2].
	@param default_val: The value that deltah(i,j) will have if no value
		is listed for that grid cell
	@param timestep_s: Length of the current coupling timestep [s]
	@param surface_temp OUT: Resulting surface temperature to use as the Dirichlet B.C.
	*/
	PetscErrorCode construct_surface_temp(
		pism::IceModelVec2S &deltah,			// IN: Input from Glint2
		double default_val,
		double timestep_s,		// Length of this coupling interval [s]
		pism::IceModelVec2S &surface_temp);

};

}}
