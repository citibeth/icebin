#pragma once

// --------------------------------
// PISM Includes... want to be included first
#include <petsc.h>
#include "IceGrid.hh"
#include "iceModelVec.hh"
// --------------------------------

namespace glint2{
namespace gpism{

/** Encapsulates mass and enthalpy together.  Used to tabulate total
enthalpy of a bunch of advected H2O based on its mass and specific
enthalpy.  This allows us to have only one C++ variable per advected
quanity, instead of two. */
struct MassEnthVec2S {
	pism::IceModelVec2S mass;
	pism::IceModelVec2S enth;

	PetscErrorCode create(pism::IceGrid &my_grid, const std::string &my_name,
		pism::IceModelVecKind ghostedp, int width = 1);


	PetscErrorCode set_attrs(
		const std::string &my_pism_intent,
		const std::string &my_long_name,
		const std::string &my_units);

#if 0
	/** @param mass [kg m-2]
	@param specific_enthalpy [J kg-1] */
	void add_mass(int i, int j, double mass, double specific_enthalpy) {
		mass(i,j) += mass;
		enthalpy(i,j) += mass * specific_enthalpy;
	}

	void clear() {
		mass.set(0);
		enth.set(0);
	}
#endif
};

struct VecWithFlags {
	pism::IceModelVec2S &vec;
	int flags;

	VecWithFlags(pism::IceModelVec2S &_vec, int _flags) :
		vec(_vec), flags(_flags) {}
};

class MassEnergyBudget {
public:
	// ============================================================
	// Total State

	// ------------ Enthalpy State
	MassEnthVec2S total;			// Total mass [kg m-2] and enthalpy [J m-2] of the ice sheet

	// =============================================================
	// Cumulative Fluxes

	// ======================= Variables to accumulate PISM output
	// These are accumulated as [kg m-2] or [J m-2]
	// Sign: positive always means it ADDS to the ice sheet, negative REMOVES from the ice sheet.
	// They are accumulated for the life of the simulation, and never zeroed.
	// Other instances of MassEnergyBudget are used to compute differences.

	// ----------- Heat generation of flows [vertical]
	pism::IceModelVec2S basal_frictional_heating;	//!< Total amount of basal friction heating [J/m^2]
	pism::IceModelVec2S strain_heating;	//!< Total amount of strain heating [J/m^2]
	pism::IceModelVec2S geothermal_flux;	//!< Total amount of geothermal energy [J/m^2]
	pism::IceModelVec2S upward_geothermal_flux;	//!< Total amount of geothermal energy [J/m^2]


	// ----------- Mass advection, with accompanying enthalpy change
	// The enthalpy reported for these mass fluxes are the enthalpies
	// AS REPORTED TO GLINT2!  That is not necessarily the same as the enthalpy
	// that PISM sees internally.
	MassEnthVec2S calving;			//!< Equal to IceModel::discharge_flux_2D_cumulative
	MassEnthVec2S basal_runoff;		//!< Enthalpy here is predictable, since runoff is 0C 100% water fraction.
	MassEnthVec2S surface_mass_balance;		//!< accumulation / ablation, as provided by Glint2
	MassEnthVec2S melt_grounded;		//!< basal melt (grounded) (from summing meltrate_grounded)
	MassEnthVec2S melt_floating;		//!< sub-shelf melt (from summing meltrate_floating)

	// ----------- Mass advection WITHIN the ice sheet
	MassEnthVec2S internal_advection;
//	MassEnthVec2S divQ_SIA;
//	MassEnthVec2S divQ_SSA;		//!< flux divergence

	// ======================= Balance the Budget
	// At each step, we set epsilon as follows:
	// total - (sum of fluxes) + epsilon = 0
	// ==> epsilon = (sum of fluxes) - total
	MassEnthVec2S epsilon;

	// ======================== Different sets (flags)
	static const int MASS = 1;
	static const int ENTH = 2;

	static const int TOTAL = 4;		// To be differenced at the end.
	static const int DELTA = 8;
	static const int EPSILON = 16;		// To be differenced at the end.

	// ======================== Summary of above variables
	// This makes it easy to difference two MassEnergyBudget instances.
	std::vector<VecWithFlags> all_vecs;

// =====================================================================

protected:
	void add_enth(pism::IceModelVec2S &vec, int flags)
		{ all_vecs.push_back(VecWithFlags(vec, ENTH | flags)); }

	void add_massenth(MassEnthVec2S &massenth, int flags) {
		all_vecs.push_back(VecWithFlags(massenth.mass, MASS | flags));
		all_vecs.push_back(VecWithFlags(massenth.enth, ENTH | flags));
	}

public:


	MassEnergyBudget();

	PetscErrorCode create(pism::IceGrid &grid, std::string const &prefix,
		pism::IceModelVecKind ghostedp, unsigned int width = 1);

	PetscErrorCode set_epsilon(pism::IceGrid &grid);
};

}}