#include <glint2/pism/MassEnergyBudget.hpp>

namespace glint2{
namespace gpism{


PetscErrorCode MassEnthVec2S::create(pism::IceGrid &my_grid, const std::string &my_name,
	pism::IceModelVecKind ghostedp, int width)
{
	PetscErrorCode ierr;

	ierr = mass.create(my_grid, my_name + ".mass", ghostedp, width); CHKERRQ(ierr);
	ierr = enth.create(my_grid, my_name + ".enth", ghostedp, width); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode MassEnthVec2S::set_attrs(
	const std::string &my_pism_intent,
	const std::string &my_long_name,
	const std::string &my_units)
{
	PetscErrorCode ierr;

	ierr = mass.set_attrs(my_pism_intent, my_long_name + " (mass portion)",
		"kg " + my_units, "", 0); CHKERRQ(ierr);
	ierr = enth.set_attrs(my_long_name + " (enthalpy portion)",
		"J " + my_units, "", 0); CHKERRQ(ierr);

	return 0;
}

MassEnergyBudget::MassEnergyBudget()
{
	// Energy deltas
	add_enth(basal_frictional_heating, DELTA);
	add_enth(basal_frictional_heating, DELTA);
	add_enth(strain_heating, DELTA);	//!< Total amount of strain heating [J/m^2]
	add_enth(geothermal_flux, DELTA);	//!< Total amount of geothermal energy [J/m^2]
	add_enth(upward_geothermal_flux, DELTA);	//!< Total amount of geothermal energy [J/m^2]

	// Mass and enthalpy deltas
	add_massenth(calving, DELTA);
	add_massenth(basal_runoff, DELTA);
	add_massenth(surface_mass_balance, DELTA);
	add_massenth(melt_grounded, DELTA);
	add_massenth(melt_floating, DELTA);

	// ----------- Mass advection WITHIN the ice sheet
	add_massenth(internal_advection, DELTA);

	add_massenth(epsilon, DELTA);

}

PetscErrorCode MassEnergyBudget::create(pism::IceGrid &grid, std::string const &prefix,
	pism::IceModelVecKind ghostedp, unsigned int width)
{
	PetscErrorCode ierr;

	// ----------- Mass and Enthalpy State of the Ice Sheet
	ierr = total.create(grid, prefix+"total",
		ghostedp, width); CHKERRQ(ierr);
	ierr = total.set_attrs("diagnostic",
		"State of the ice sheet (NOT a difference between timetseps)",
		"m-2 s-1"); CHKERRQ(ierr);

	// ----------- Heat generation of flows [vertical]
	// Postive means heat is flowing INTO the ice sheet.
	ierr = basal_frictional_heating.create(grid, prefix+"basal_frictional_heating",
		ghostedp, width); CHKERRQ(ierr);
	ierr = basal_frictional_heating.set_attrs("internal",
		"Basal frictional heating",
		"W m-2", ""); CHKERRQ(ierr);

	ierr = strain_heating.create(grid, prefix+"strain_heating",
		ghostedp, width); CHKERRQ(ierr);
	ierr = strain_heating.set_attrs("internal",
		"Strain heating",
		"W m-2", ""); CHKERRQ(ierr);

	ierr = geothermal_flux.create(grid, prefix+"geothermal_flux",
		ghostedp, width); CHKERRQ(ierr);
	ierr = geothermal_flux.set_attrs("internal",
		"Geothermal energy through (compare to upward_geothermal_flux?)",
		"W m-2", ""); CHKERRQ(ierr);

	ierr = upward_geothermal_flux.create(grid, prefix+"upward_geothermal_flux",
		ghostedp, width); CHKERRQ(ierr);
	ierr = upward_geothermal_flux.set_attrs("internal",
		"Geothermal energy through (compare to geothermal_flux?)",
		"W m-2", ""); CHKERRQ(ierr);

	// ----------- Mass advection, with accompanying enthalpy change
	// Postive means mass/enthalpy is flowing INTO the ice sheet.
	ierr = calving.create(grid, prefix+"calving",
		ghostedp, width); CHKERRQ(ierr);
	ierr = calving.set_attrs("diagnostic",
		"Mass/Enthalpy gain from calving.  Should be negative.",
		"m-2 s-1"); CHKERRQ(ierr);

	ierr = basal_runoff.create(grid, prefix+"basal_runoff",
		ghostedp, width); CHKERRQ(ierr);
	ierr = basal_runoff.set_attrs("diagnostic",
		"Runoff from base, should be negative.  Enthalpy portion is predictable, since runoff is 0C 100% water fraction.",
		"m-2 s-1"); CHKERRQ(ierr);

	ierr = surface_mass_balance.create(grid, prefix+"surface_mass_balance",
		ghostedp, width); CHKERRQ(ierr);
	ierr = surface_mass_balance.set_attrs("diagnostic",
		"surface_mass_balance",
		"m-2 s-1"); CHKERRQ(ierr);

	ierr = melt_grounded.create(grid, prefix+"melt_grounded",
		ghostedp, width); CHKERRQ(ierr);
	ierr = melt_grounded.set_attrs("diagnostic",
		"Basal melting of grounded ice (negative)",
		"m-2 s-1"); CHKERRQ(ierr);

	ierr = melt_floating.create(grid, prefix+"melt_floating",
		ghostedp, width); CHKERRQ(ierr);
	ierr = melt_floating.set_attrs("diagnostic",
		"Sub-shelf melting (negative)",
		"m-2 s-1"); CHKERRQ(ierr);

	// ----------- Advection WITHIN the ice sheet
	ierr = internal_advection.create(grid, prefix+"internal_advection",
		ghostedp, width); CHKERRQ(ierr);
	ierr = internal_advection.set_attrs("diagnostic",
		"Advection within the ice sheet",
		"m-2 s-1"); CHKERRQ(ierr);

	// ----------- Balance the Budget
	ierr = epsilon.create(grid, prefix+"epsilon",
		ghostedp, width); CHKERRQ(ierr);
	ierr = epsilon.set_attrs("diagnostic",
		"Unaccounted-for changes",
		"m-2 s-1"); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode MassEnergyBudget::set_epsilon(pism::IceGrid &grid)
{
	PetscErrorCode ierr;

	// ==> epsilon = (sum of fluxes) - total

	// -------- Mass
	ierr = epsilon.mass.begin_access(); CHKERRQ(ierr);
	ierr = total.mass.begin_access(); CHKERRQ(ierr);
	for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
	for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
		epsilon.mass(i,j) = -total.mass(i,j);
	}}
	ierr = total.mass.end_access(); CHKERRQ(ierr);

	for (auto &ii : all_vecs) {
		if (!(ii.flags & (DELTA | MASS))) continue;

		ierr = ii.vec.begin_access(); CHKERRQ(ierr);
		for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
		for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
			epsilon.mass(i,j) += ii.vec(i,j);
		}}
		ierr = ii.vec.end_access(); CHKERRQ(ierr);
	}
	ierr = epsilon.mass.end_access(); CHKERRQ(ierr);

	// -------- Energy
	ierr = epsilon.enth.begin_access(); CHKERRQ(ierr);
	ierr = total.enth.begin_access(); CHKERRQ(ierr);
	for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
	for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
		epsilon.enth(i,j) = -total.enth(i,j);
	}}
	ierr = total.enth.end_access(); CHKERRQ(ierr);

	for (auto &ii : all_vecs) {
		if (!(ii.flags & (DELTA | ENTH))) continue;

		for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
		for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
			epsilon.enth(i,j) += ii.vec(i,j);
		}}
		ierr = ii.vec.end_access(); CHKERRQ(ierr);
	}
	ierr = epsilon.enth.end_access(); CHKERRQ(ierr);
	return 0;
}

}}
