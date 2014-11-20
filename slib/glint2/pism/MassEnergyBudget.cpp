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
	const std::string &my_units,
	const std::string &my_standard_name)
{
	PetscErrorCode ierr;

	ierr = mass.set_attrs(my_pism_intent, my_long_name + " (mass portion)",
		"kg " + my_units, my_standard_name + ".mass", 0); CHKERRQ(ierr);
	ierr = enth.set_attrs(my_pism_intent, my_long_name + " (enthalpy portion)",
		"J " + my_units, my_standard_name + ".enth", 0); CHKERRQ(ierr);

	return 0;
}

MassEnergyBudget::MassEnergyBudget()
{
}

PetscErrorCode MassEnergyBudget::create(pism::IceGrid &grid, std::string const &prefix,
	pism::IceModelVecKind ghostedp, unsigned int width)
{
	PetscErrorCode ierr;

	printf("BEGIN MassEnergyBudget::create()\n");

	// ----------- Mass and Enthalpy State of the Ice Sheet
	ierr = total.create(grid, prefix+"total",
		ghostedp, width); CHKERRQ(ierr);
	ierr = total.set_attrs("diagnostic",
		"State of the ice sheet (NOT a difference between timetseps)",
		"m-2", "total"); CHKERRQ(ierr);
	add_massenth(total, TOTAL, "", "");

	// ----------- Heat generation of flows [vertical]
	// Postive means heat is flowing INTO the ice sheet.
	ierr = basal_frictional_heating.create(grid, prefix+"basal_frictional_heating",
		ghostedp, width); CHKERRQ(ierr);
	ierr = basal_frictional_heating.set_attrs("internal",
		"Basal frictional heating",
		"W m-2", ""); CHKERRQ(ierr);
	add_enth(basal_frictional_heating, DELTA, "basal_frictional_heating");

	ierr = strain_heating.create(grid, prefix+"strain_heating",
		ghostedp, width); CHKERRQ(ierr);
	ierr = strain_heating.set_attrs("internal",
		"Strain heating",
		"W m-2", ""); CHKERRQ(ierr);
	add_enth(strain_heating, DELTA, "strain_heating");	//!< Total amount of strain heating [J/m^2]

	ierr = geothermal_flux.create(grid, prefix+"geothermal_flux",
		ghostedp, width); CHKERRQ(ierr);
	ierr = geothermal_flux.set_attrs("internal",
		"Geothermal energy through (compare to upward_geothermal_flux?)",
		"W m-2", ""); CHKERRQ(ierr);
	add_enth(geothermal_flux, DELTA, "geothermal_flux");	//!< Total amount of geothermal energy [J/m^2]

	ierr = upward_geothermal_flux.create(grid, prefix+"upward_geothermal_flux",
		ghostedp, width); CHKERRQ(ierr);
	ierr = upward_geothermal_flux.set_attrs("internal",
		"Geothermal energy through (compare to geothermal_flux?)",
		"W m-2", ""); CHKERRQ(ierr);
	add_enth(upward_geothermal_flux, DELTA, "upward_geothermal_flux");	//!< Total amount of geothermal energy [J/m^2]

	// ----------- Mass advection, with accompanying enthalpy change
	// Postive means mass/enthalpy is flowing INTO the ice sheet.
	std::string name;

	ierr = calving.create(grid, prefix+"calving",
		ghostedp, width); CHKERRQ(ierr);
	ierr = calving.set_attrs("diagnostic",
		"Mass/Enthalpy gain from calving.  Should be negative.",
		"m-2 s-1", "calving"); CHKERRQ(ierr);
	add_massenth(calving, DELTA, "calving.mass", "calving.enth");

#if 0
	ierr = basal_runoff.create(grid, prefix+"basal_runoff",
		ghostedp, width); CHKERRQ(ierr);
	ierr = basal_runoff.set_attrs("diagnostic",
		"Runoff from base, should be negative.  Enthalpy portion is predictable, since runoff is 0C 100% water fraction.",
		"m-2 s-1"); CHKERRQ(ierr);
	add_massenth(basal_runoff, DELTA, "basal_runoff.mass", "basal_runoff.enth");
#endif

	ierr = surface_mass_balance.create(grid, prefix+"surface_mass_balance",
		ghostedp, width); CHKERRQ(ierr);
	ierr = surface_mass_balance.set_attrs("diagnostic",
		"surface_mass_balance",
		"m-2 s-1", "surface_mass_balance"); CHKERRQ(ierr);
	add_massenth(surface_mass_balance, DELTA, "surface_mass_balance.mass", "surface_mass_balance.enth");

	ierr = pism_smb.create(grid, prefix+"pism_smb",
		ghostedp, width); CHKERRQ(ierr);
	ierr = pism_smb.set_attrs("diagnostic",
		"pism_smb",
		"kg m-2 s-1", "pism_smb"); CHKERRQ(ierr);
	// No DELTA< does not participate in epsilon computation
	add_mass(pism_smb, 0, "", "");

	ierr = href_to_h.create(grid, prefix+"href_to_h",
		ghostedp, width); CHKERRQ(ierr);
	ierr = href_to_h.set_attrs("diagnostic",
		"href_to_h",
		"kg m-2 s-1", "href_to_h"); CHKERRQ(ierr);
	add_mass(href_to_h, 0, "");

	ierr = nonneg_rule.create(grid, prefix+"nonneg_rule",
		ghostedp, width); CHKERRQ(ierr);
	ierr = nonneg_rule.set_attrs("diagnostic",
		"nonneg_rule",
		"kg m-2 s-1", "nonneg_rule"); CHKERRQ(ierr);
	add_mass(nonneg_rule, 0, "");



	ierr = melt_grounded.create(grid, prefix+"melt_grounded",
		ghostedp, width); CHKERRQ(ierr);
	ierr = melt_grounded.set_attrs("diagnostic",
		"Basal melting of grounded ice (negative)",
		"m-2 s-1", "melt_grounded"); CHKERRQ(ierr);
	add_massenth(melt_grounded, DELTA, "melt_grounded.mass", "melt_grounded.enth");

	ierr = melt_floating.create(grid, prefix+"melt_floating",
		ghostedp, width); CHKERRQ(ierr);
	ierr = melt_floating.set_attrs("diagnostic",
		"Sub-shelf melting (negative)",
		"m-2 s-1", "melt_floating"); CHKERRQ(ierr);
	add_massenth(melt_floating, DELTA, "melt_floating.mass", "melt_floating.enth");

	// ----------- Advection WITHIN the ice sheet
	ierr = internal_advection.create(grid, prefix+"internal_advection",
		ghostedp, width); CHKERRQ(ierr);
	ierr = internal_advection.set_attrs("diagnostic",
		"Advection within the ice sheet",
		"m-2 s-1", "internal_advection"); CHKERRQ(ierr);
	add_massenth(internal_advection, DELTA, "internal_advection.mass", "internal_advection.enth");

	// ----------- Balance the Budget
	ierr = epsilon.create(grid, prefix+"epsilon",
		ghostedp, width); CHKERRQ(ierr);
	ierr = epsilon.set_attrs("diagnostic",
		"Unaccounted-for changes",
		"m-2 s-1", "epsilon"); CHKERRQ(ierr);
	add_massenth(epsilon, EPSILON, "epsilon.mass", "epsilon.enth");

	return 0;
}

PetscErrorCode MassEnergyBudget::set_epsilon(pism::IceGrid &grid)
{
	PetscErrorCode ierr;

	// ==> epsilon = (sum of fluxes) - total
	printf("BEGIN MassEnergyBudget::set_epsilon()\n");

	// -------- Mass
	ierr = epsilon.mass.begin_access(); CHKERRQ(ierr);
	ierr = total.mass.begin_access(); CHKERRQ(ierr);
	for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
	for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
		epsilon.mass(i,j) = total.mass(i,j);
	}}
	ierr = total.mass.end_access(); CHKERRQ(ierr);

	for (auto &ii : all_vecs) {
		if ((ii.flags & (DELTA | MASS)) != (DELTA | MASS)) continue;

		printf("epsilon.mass: %s\n", ii.vec.name().c_str());

		ierr = ii.vec.begin_access(); CHKERRQ(ierr);
		for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
		for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
			epsilon.mass(i,j) -= ii.vec(i,j);
		}}
		ierr = ii.vec.end_access(); CHKERRQ(ierr);
	}
	ierr = epsilon.mass.end_access(); CHKERRQ(ierr);

	// -------- Energy
	ierr = epsilon.enth.begin_access(); CHKERRQ(ierr);
	ierr = total.enth.begin_access(); CHKERRQ(ierr);
	for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
	for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
		epsilon.enth(i,j) = total.enth(i,j);
	}}
	ierr = total.enth.end_access(); CHKERRQ(ierr);

#if 1
	for (auto &ii : all_vecs) {
		if ((ii.flags & (DELTA | ENTH)) != (DELTA | ENTH)) continue;

		printf("epsilon.energy: %s\n", ii.vec.name().c_str());

		ierr = ii.vec.begin_access(); CHKERRQ(ierr);
		for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
		for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
			epsilon.enth(i,j) -= ii.vec(i,j);
		}}
		ierr = ii.vec.end_access(); CHKERRQ(ierr);
	}
#endif
	ierr = epsilon.enth.end_access(); CHKERRQ(ierr);

	printf("END MassEnergyBudget::set_epsilon()\n");
	return 0;
}

}}
