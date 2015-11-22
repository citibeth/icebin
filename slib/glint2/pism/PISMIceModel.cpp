#include <cfloat>
#include <iostream>
#include <glint2/pism/PISMIceModel.hpp>
#include <glint2/pism/GLINT2EnthalpyConverter.hpp>
#include "bedrockThermalUnit.hh"

using namespace pism;

namespace glint2 {
namespace gpism {

// /** Sum a 3-D vector in the Z direction to create a 2-D vector.
// 
// <p>Note that this sums up all the values in a column, including ones
// above the ice. This may or may not be what you need. Also, take a look
// at IceModel::compute_ice_enthalpy(PetscScalar &result) in iMreport.cc.</p>
// 
// <p>As for the difference between IceModelVec2 and IceModelVec2S, the
// former can store fields with more than 1 "degree of freedom" per grid
// point (such as 2D fields on the "staggered" grid, with the first
// degree of freedom corresponding to the i-offset and second to
// j-offset).</p>
// 
// <p>IceModelVec2S is just IceModelVec2 with "dof == 1", and
// IceModelVec2V is IceModelVec2 with "dof == 2". (Plus some extra
// methods, of course.)</p>
// 
// <p>Either one of IceModelVec2 and IceModelVec2S would work in this
// case.</p>
// 
// @see https://github.com/pism/pism/issues/229 */
// static PetscErrorCode sum_columns(IceModelVec3 &input, IceModelVec2S &output)
// {
//   PetscScalar *column = NULL;
//   IceGrid &grid = *input.get_grid();
//   int ierr = 0;
// 
//   ierr = input.begin_access(); CHKERRQ(ierr);
//   ierr = output.begin_access(); CHKERRQ(ierr);
//   for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
//     for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
//       ierr = input.getInternalColumn(i, j, &column); CHKERRQ(ierr);
// 
//       output(i,j) = 0.0;
//       for (unsigned int k = 0; k < grid.Mz; ++k)
//         output(i,j) += column[k];
//     }
//   }
//   ierr = output.end_access(); CHKERRQ(ierr);
//   ierr = input.end_access(); CHKERRQ(ierr);
// 
//   return 0;
// }
// 

// ================================

PISMIceModel::PISMIceModel(pism::IceGrid &g, pism::Config &config, pism::Config &overrides, PISMIceModel::Params const &_params) :
	pism::IceModel(g, config, overrides),
	params(_params)
{}

PISMIceModel::~PISMIceModel() {} // must be virtual merely because some members are virtual


//! \brief Decide which enthalpy converter to use.
PetscErrorCode PISMIceModel::allocate_enthalpy_converter() {
	PetscErrorCode ierr;

	if (EC != NULL)
		return 0;

	EC = new GLINT2EnthalpyConverter(config);

	if (getVerbosityLevel() > 3) {
		PetscViewer viewer;
		ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer); CHKERRQ(ierr);
		ierr = EC->viewConstants(viewer); CHKERRQ(ierr);
	}

	return 0;
}

PetscErrorCode PISMIceModel::allocate_subglacial_hydrology()
{
	printf("BEGIN PISMIceModel::allocate_subglacial_hydrology()\n");
//	printf("subglacial_hydrology = %p %p\n", subglacial_hydrology, pism::IceModel::subglacial_hydrology);
	if (pism::IceModel::subglacial_hydrology != NULL) return 0; // indicates it has already been allocated
	subglacial_hydrology = new NullTransportHydrology(grid, config);
	printf("END PISMIceModel::allocate_subglacial_hydrology()\n");
	return 0;
}



PetscErrorCode PISMIceModel::allocate_couplers()
{
	printf("BEGIN PISMIceModel::allocate_couplers()\n");
	PetscErrorCode ierr;
	// Initialize boundary models:
	PAFactory pa(grid, config);
	PSFactory ps(grid, config);
	POFactory po(grid, config);
	pism::AtmosphereModel *atmosphere;

	ierr = PetscOptionsBegin(grid.com, "", "Options choosing PISM boundary models", ""); CHKERRQ(ierr);

	// GLINT2-modified version
	if (surface == NULL) {
		surface = new PSConstantGLINT2(grid, config);
		external_surface_model = false;

		pa.create(atmosphere);
		surface->attach_atmosphere_model(atmosphere);
	}

	if (ocean == NULL) {
		po.create(ocean);
		external_ocean_model = false;
	}
	ierr = PetscOptionsEnd(); CHKERRQ(ierr);

	printf("END PISMIceModel::allocate_couplers()\n");
	return 0;
}


PetscErrorCode PISMIceModel::createVecs()
{
	super::createVecs();

	PetscErrorCode ierr;

	printf("BEGIN PISMIceModel::createVecs()\n");
	ierr = base.create(grid, "", WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = cur.create(grid, "", WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = rate.create(grid, "", WITHOUT_GHOSTS); CHKERRQ(ierr);
	printf("END PISMIceModel::createVecs()\n");

	ierr = mask.create(grid, "mask", pism::WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = mask.set_attrs("diagnostic", "Land surface types", "", "mask");

	ierr = M1.create(grid, "M1", pism::WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = M1.set_attrs("diagnostic", "Mass of top layer", "kg m-2", "M1");
	ierr = M2.create(grid, "M2", pism::WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = M2.set_attrs("diagnostic", "Mass of second-top layer", "kg m-2", "M2");

	ierr = H1.create(grid, "H1", pism::WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = H1.set_attrs("diagnostic", "Enthalpy of top layer", "J m-2", "H1");
	ierr = H2.create(grid, "H2", pism::WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = H2.set_attrs("diagnostic", "Enthalpy of second-top layer", "J m-2", "H2");

	ierr = V1.create(grid, "V1", pism::WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = V1.set_attrs("diagnostic", "Volume of top layer", "m^3 m-2", "V1");
	ierr = V2.create(grid, "V2", pism::WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = V2.set_attrs("diagnostic", "Volume of second-top layer", "m^3 m-2", "V2");

	std::cout << "PISMIceModel Conservation Formulas:" << std::endl;
	cur.print_formulas(std::cout);

	return 0;
}

PetscErrorCode PISMIceModel::massContPreHook()
{
#if 0	// Avoid warnings until we put something in this method
	PetscErrorCode ierr;

	// Enthalpy and mass continuity are stepped with different timesteps.
	// Fish out the timestep relevant to US.
	const double my_t0 = grid.time->current();
	const double my_dt = this->dt;
#endif

	return 0;
}


PetscErrorCode PISMIceModel::massContPostHook()
{
#if 0	// Avoid warnings until we put something in this method
	PetscErrorCode ierr;

	// Enthalpy and mass continuity are stepped with different timesteps.
	// Fish out the timestep relevant to US.
	const double my_t0 = grid.time->current();
	const double my_dt = this->dt;
#endif


	return 0;
}


PetscErrorCode PISMIceModel::energyStep()
{
	PetscErrorCode ierr;

	printf("BEGIN PISMIceModel::energyStep(t=%f, dt=%f)\n", t_TempAge, dt_TempAge);

	// Enthalpy and mass continuity are stepped with different timesteps.
	// Fish out the timestep relevant to US.
	// const double my_t0 = t_TempAge;			// Time at beginning of timestep
	const double my_dt = dt_TempAge;

	// =========== BEFORE Energy Step

	// =========== The Energy Step Itself
	ierr = super::energyStep(); CHKERRQ(ierr);

	// =========== AFTER Energy Step

	// We need to integrate over strain_heating and geothermal_flux, which
	// are given in PISM as rates.


	// --------- Upward Geothermal Flux
	// Use actual geothermal flux, not the long-term average..
	// See: file:///Users/rpfische/git/pism/build/doc/browser/html/classPISMBedThermalUnit.html#details
	{
		IceModelVec2S &upward_geothermal_flux(vWork2d[0]);
		ierr = btu->get_upward_geothermal_flux(upward_geothermal_flux); CHKERRQ(ierr);
		ierr = cur.upward_geothermal_flux.add(my_dt, upward_geothermal_flux); CHKERRQ(ierr);
	}

	// ----------- Geothermal Flux
	ierr = cur.geothermal_flux.add(my_dt, geothermal_flux); CHKERRQ(ierr);

	// ---------- Basal Frictional Heating (see iMenthalpy.cc l. 220)
	IceModelVec2S *Rb = NULL;
	ierr = stress_balance->get_basal_frictional_heating(Rb); CHKERRQ(ierr);
	cur.basal_frictional_heating.add(my_dt, *Rb);

	// NOTE: strain_heating is inf at the coastlines.
	// See: https://github.com/pism/pism/issues/292
	// ------------ Volumetric Strain Heating
	// strain_heating_sum += my_dt * sum_columns(strainheating3p)
	IceModelVec3 *strain_heating3p;
	stress_balance->get_volumetric_strain_heating(strain_heating3p);
	// cur.strain_heating = cur.strain_heating * 1.0 + my_dt * sum_columns(strain_heating3p)
	ierr = strain_heating3p->sumColumns(cur.strain_heating, 1.0, my_dt); CHKERRQ(ierr);

	printf("END PISMIceModel::energyStep(time=%f)\n", t_TempAge);
	return 0;
}

PetscErrorCode PISMIceModel::massContExplicitStep() {
	PetscErrorCode ierr;

	printf("BEGIN PISMIceModel::MassContExplicitStep(dt = %g)\n", dt);

	_ice_density = config.get("ice_density");
	_meter_per_s_to_kg_per_m2 = dt * _ice_density;


	// =========== The Mass Continuity Step Itself
	// This will call through to accumulateFluxes_massContExplicitStep()
	// in the inner loop
	ierr = Enth3.begin_access(); CHKERRQ(ierr);
	ierr = cur.pism_smb.begin_access(); CHKERRQ(ierr);
	ierr = cur.melt_grounded.begin_access(); CHKERRQ(ierr);
	ierr = cur.melt_floating.begin_access(); CHKERRQ(ierr);
	ierr = cur.internal_advection.begin_access(); CHKERRQ(ierr);
	ierr = cur.href_to_h.begin_access(); CHKERRQ(ierr);
	ierr = cur.nonneg_rule.begin_access(); CHKERRQ(ierr);

	ierr = super::massContExplicitStep(); CHKERRQ(ierr);

	ierr = cur.nonneg_rule.end_access(); CHKERRQ(ierr);
	ierr = cur.href_to_h.end_access(); CHKERRQ(ierr);
	ierr = cur.internal_advection.end_access(); CHKERRQ(ierr);
	ierr = cur.melt_floating.end_access(); CHKERRQ(ierr);
	ierr = cur.melt_grounded.end_access(); CHKERRQ(ierr);
	ierr = cur.pism_smb.end_access(); CHKERRQ(ierr);
	ierr = Enth3.end_access(); CHKERRQ(ierr);



	// =========== AFTER the Mass Continuity Step

	// ----------- SMB: Pass inputs through to outputs.
	// They are needed to participate in mass/energy budget
	PSConstantGLINT2 *surface = ps_constant_glint2();
	ierr = surface->glint2_massxfer_rate.begin_access(); CHKERRQ(ierr);
	ierr = surface->glint2_enthxfer_rate.begin_access(); CHKERRQ(ierr);
	ierr = surface->glint2_deltah.begin_access(); CHKERRQ(ierr);
	ierr = cur.glint2_smb.begin_access(); CHKERRQ(ierr);
	ierr = cur.glint2_deltah.begin_access(); CHKERRQ(ierr);
	for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
	for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
		cur.glint2_smb.mass(i,j) += dt * surface->glint2_massxfer_rate(i,j);
		cur.glint2_smb.enth(i,j) += dt * surface->glint2_enthxfer_rate(i,j);
		cur.glint2_deltah(i,j) += surface->glint2_deltah(i,j);
	}}
	ierr = surface->glint2_massxfer_rate.end_access(); CHKERRQ(ierr);
	ierr = surface->glint2_enthxfer_rate.end_access(); CHKERRQ(ierr);
	ierr = surface->glint2_deltah.end_access(); CHKERRQ(ierr);
	ierr = cur.glint2_smb.end_access(); CHKERRQ(ierr);
	ierr = cur.glint2_deltah.end_access(); CHKERRQ(ierr);

	printf("END PISMIceModel::MassContExplicitStep()\n");
	return 0;
}



/** This is called IMMEDIATELY after ice is gained/lost in
iMgeometry.cc (massContExplicitStep()).  Here we can record the same
values that PISM saw when moving ice around. */
PetscErrorCode PISMIceModel::accumulateFluxes_massContExplicitStep(
	int i, int j,
	double surface_mass_balance,		// [m s-1] ice equivalent
	double meltrate_grounded,			// [m s-1] ice equivalent
	double meltrate_floating,			// [m s-1] ice equivalent
	double divQ_SIA,					// [m s-1] ice equivalent
	double divQ_SSA,					// [m s-1] ice equivalent
	double Href_to_H_flux,				// [m] ice equivalent
	double nonneg_rule_flux)			// [m] ice equivalent
{
	PetscErrorCode ierr;

	// -------------- Melting
	double p_basal = EC->getPressureFromDepth(ice_thickness(i,j));
	double T = EC->getMeltingTemp(p_basal);
	double specific_enth_basal;
	ierr = EC->getEnthPermissive(T, 1.0, p_basal, specific_enth_basal); CHKERRQ(ierr);
	double mass;

	// ------- Melting at base of ice sheet
	mass = -meltrate_grounded * _meter_per_s_to_kg_per_m2;
	cur.melt_grounded.mass(i,j) += mass;
	cur.melt_grounded.enth(i,j) += mass * specific_enth_basal;

	// ------- Melting under ice shelf
	mass = -meltrate_floating * _meter_per_s_to_kg_per_m2;
	cur.melt_floating.mass(i,j) += mass;
	cur.melt_floating.enth(i,j) += mass * specific_enth_basal;


	// -------------- internal_advection
	const int ks = grid.kBelowHeight(ice_thickness(i,j));
	double *Enth;
	ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
	double specific_enth_top = Enth[ks];		// Approximate, we will use the enthalpy of the top layer...

	mass = -(divQ_SIA + divQ_SSA) * _meter_per_s_to_kg_per_m2;

	cur.internal_advection.mass(i,j) += mass;
	cur.internal_advection.enth(i,j) += mass * specific_enth_top;


	// -------------- Get the easy veriables out of the way...
	mass = surface_mass_balance * _meter_per_s_to_kg_per_m2;
	cur.pism_smb.mass(i,j) += mass;
	cur.pism_smb.enth(i,j) += mass * specific_enth_top;
	cur.nonneg_rule(i,j) -= nonneg_rule_flux * _ice_density;
	cur.href_to_h(i,j) += Href_to_H_flux * _ice_density;


//	printf("END PISMIceModel::accumulateFluxes_MassContExplicitStep()\n");
	return 0;
}


/** @param t0 Time of last time we coupled. */
PetscErrorCode PISMIceModel::set_rate(double dt)
{
	PetscErrorCode ierr;

	printf("BEGIN PISMIceModel::set_rate(dt=%f)\n", dt);

	double by_dt = 1.0 / dt;

	ierr = compute_enth2(cur.total.enth, cur.total.mass); CHKERRQ(ierr);
	ierr = cur.set_epsilon(grid); CHKERRQ(ierr);

	// Compute differences, and set base = cur
	auto base_ii(base.all_vecs.begin());
	auto cur_ii(cur.all_vecs.begin());
	auto rate_ii(rate.all_vecs.begin());
	for (; base_ii != base.all_vecs.end(); ++base_ii, ++cur_ii, ++rate_ii) {
		IceModelVec2S &vbase(base_ii->vec);
		IceModelVec2S &vcur(cur_ii->vec);
		IceModelVec2S &vrate(rate_ii->vec);

		ierr = vbase.begin_access(); CHKERRQ(ierr);
		ierr = vcur.begin_access(); CHKERRQ(ierr);
		ierr = vrate.begin_access(); CHKERRQ(ierr);
		for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
		for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
			// rate = cur - base: Just for DELTA and EPISLON flagged vectors
			if (base_ii->flags & (MassEnergyBudget::DELTA | MassEnergyBudget::EPSILON)) {
				vrate(i,j) = (vcur(i,j) - vbase(i,j)) * by_dt;
			} else {
				// Or else just copy the to rate
				vrate(i,j) = vcur(i,j);
			}
		}}
		ierr = vrate.end_access(); CHKERRQ(ierr);
		ierr = vcur.end_access(); CHKERRQ(ierr);
		ierr = vbase.end_access(); CHKERRQ(ierr);

	}

	printf("END PISMIceModel::set_rate()\n");
	return 0;

}

PetscErrorCode PISMIceModel::reset_rate()
{
	PetscErrorCode ierr;


	// Compute differences, and set base = cur
	auto base_ii(base.all_vecs.begin());
	auto cur_ii(cur.all_vecs.begin());
	for (; base_ii != base.all_vecs.end(); ++base_ii, ++cur_ii) {
		IceModelVec2S &vbase(base_ii->vec);
		IceModelVec2S &vcur(cur_ii->vec);

		// This cannot go in the loop above with PETSc because
		// vbase is needed on the RHS of the equations above.
		ierr = vbase.begin_access(); CHKERRQ(ierr);
		ierr = vcur.begin_access(); CHKERRQ(ierr);
		for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
		for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
			// base = cur: For ALL vectors
			vbase(i,j) = vcur(i,j);
		}}
		ierr = vcur.end_access(); CHKERRQ(ierr);
		ierr = vbase.end_access(); CHKERRQ(ierr);
	}


#if 0
ierr = rate.geothermal_flux.begin_access(); CHKERRQ(ierr);
printf("GG rate.geothermal_flux(%d, %d) = %f (%p)\n", grid.xs, grid.xs, rate.geothermal_flux(grid.xs, grid.xs), &rate.geothermal_flux(grid.xs, grid.xs));
ierr = rate.geothermal_flux.end_access(); CHKERRQ(ierr);

ierr = cur.geothermal_flux.begin_access(); CHKERRQ(ierr);
printf("GG cur.geothermal_flux(%d, %d) = %f (%p)\n", grid.xs, grid.xs, cur.geothermal_flux(grid.xs, grid.xs), &cur.geothermal_flux(grid.xs, grid.xs));
ierr = cur.geothermal_flux.end_access(); CHKERRQ(ierr);

ierr = base.geothermal_flux.begin_access(); CHKERRQ(ierr);
printf("GG base.geothermal_flux(%d, %d) = %f (%p)\n", grid.xs, grid.xs, base.geothermal_flux(grid.xs, grid.xs), &base.geothermal_flux(grid.xs, grid.xs));
ierr = base.geothermal_flux.end_access(); CHKERRQ(ierr);
#endif

	return 0;
}

/** @param t0 Time of last time we coupled. */
PetscErrorCode PISMIceModel::prepare_outputs(double t0)
{
	PetscErrorCode ierr;

	printf("BEGIN PISMIceModel::prepare_outputs()\n");

	// ------ Difference between now and the last time we were called
	double t1 = enthalpy_t();	// Current time of the enthalpy portion of ice model.
	ierr = set_rate(t1 - t0); CHKERRQ(ierr);

	// ice_surface_enth & ice_surfac_enth_depth
	ierr = prepare_initial_outputs(); CHKERRQ(ierr);

	printf("END PISMIceModel::prepare_outputs()\n");
	return 0;
}
// -----------------------------------------------------------------------
/** This is an optimized version of the following function, optimized for two layers only:
def integrate_weights(dz, z0, z1):
	"""dz: Depth of each layer"""

	top = np.cumsum(dz)		# Top of each layer
	ix = []
	weights = []

	# Nothing to do here
	z0 = min(z0, top[-1])
	z1 = min(z1, top[-1])
	if (z0 >= z1):
		return ix,weights

	i0 = bisect.bisect_right(top, z0)
	i1 = bisect.bisect_left(top, z1)

	# First interval
	ix.append(i0)
	weights.append((min(z1, top[i0]) - z0)/dz[i0])

	if i1 > i0:
		# Middle intervals
		if i1 > i0+1:
			for i in range(i0+1,i1):
				ix.append(i)
				weights.append(1.)

		# Last interval
		ix.append(i1)
		weights.append((z1 - top[i1-1]) / dz[i1])

#	print(dz,z0,z1,ix,weights,wsum,i0,i1)
	return ix,weights
*/
static PetscErrorCode integrate_weights_2(
double dz[2],		// Thickness of each PISM layer
double z0,			// Start of integration
double z1,			// End of integration
double w[2])		// Weight of each layer
{
	w[0] = 0;
	w[1] = 0;
	double const top[2] = {dz[0], dz[0]+dz[1]};

	// Assume zero outside our range
	z0 = std::max(std::min(z0, top[1]), 0.);
	z1 = std::max(std::min(z1, top[1]), 0.);

	// Assume z0 <= z1
	if (z0 > z1) {
		fprintf(stderr, "Must have z0 (%g) <= z1 (%g)\n", z0, z1);
		return 1;
	}
	if (z0 == z1) return 0;	// All zeros now


	// ------ Figure out which grid cell z0 and z1 fall into
	// i0 = bisect.bisect_right(top, z0)
	int i0 = (z0 >= top[0] ? 1 : 0);
	int i1 = (z1 > top[0] ? 1 : 0);

	// ------- First interval
	w[i0] = (std::min(z1, top[i0]) - z0) / dz[i0];
	if (i1 > i0) {
		w[1] = (z1 - top[0]) / dz[1];
	}

	return 0;
}
// -----------------------------------------------------------------------
#if 0
static inline void set_gcm_layer(double vv[2], double ee[2], double w[2], double &V1, double &M1, double &H1, double const gcm_thk, double const ice_density)
{
	V1 = vv[0]*w[0] + vv[1]*w[1];
	M1 = V1 * ice_density;
	H1 = (ee[0]*vv[0]*w[0] + ee[1]*vv[1]*w[1]) * ice_density;

	// Extend layer to proper thickness
	if (V1 < gcm_thk) {
		double v1old = V1;
		V1 = gcm_thk;	// Extend layer to proper thickness
		double fact = V1 / v1old;
		M1 *= fact;
		H1 *= fact;
	}
}
#endif

PetscErrorCode PISMIceModel::prepare_initial_outputs()
{
	PetscErrorCode ierr;

	double ice_density = config.get("ice_density");	// [kg m-3]

	// --------- mask
	ierr = vMask.begin_access(); CHKERRQ(ierr);
	ierr = mask.begin_access(); CHKERRQ(ierr);
	for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
	for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
		mask(i,j) = vMask(i,j);
	}}
	ierr = mask.end_access(); CHKERRQ(ierr);
	ierr = vMask.end_access(); CHKERRQ(ierr);

	// --------- ice_surface_enth from Enth3
	ierr = Enth3.begin_access(); CHKERRQ(ierr);
	ierr = M1.begin_access(); CHKERRQ(ierr);
	ierr = M2.begin_access(); CHKERRQ(ierr);
	ierr = H1.begin_access(); CHKERRQ(ierr);
	ierr = H2.begin_access(); CHKERRQ(ierr);
	ierr = V1.begin_access(); CHKERRQ(ierr);
	ierr = V2.begin_access(); CHKERRQ(ierr);
	ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
	for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
	for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {

		// ----------- Obtain top two layers of PISM
		double *Enth;
		ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);

		double vv[2];	// Thickness of top two layers of PISM [m^3 m-2]
		double ee[2];	// Specific Enthalpy of top two yers of PISM [J kg-1]

		// ---- Get depth of the top two layers
		int const ks = grid.kBelowHeight(ice_thickness(i,j));
		vv[0] = ice_thickness(i,j) - grid.zlevels[ks];	// [m^3 m-2]
		ee[0] = Enth[ks];
		int const ks2 = ks - 1;
		if (ks2 >= 0) {
			vv[1] = grid.zlevels[ks] - grid.zlevels[ks2];
			ee[1] = Enth[ks2];
		} else {
			// There is no second layer
			vv[1] = 0;
			ee[1] = 0;
		}
		double E1;


		// ---------- Resample to two layers for Glint2
		const double gcm_thk = 5.;		// Thickness to make each layer we produce [m]
			// MUST have: gcm_thk*2 <= thickness of normal PISM layer (40m)
		double w[2];

		// ---- Top layer for GCM
		ierr = integrate_weights_2(vv, 0., gcm_thk, w); CHKERRQ(ierr);
		V1(i,j) = vv[0]*w[0] + vv[1]*w[1];
		M1(i,j) = V1(i,j) * ice_density;
		H1(i,j) = (ee[0]*vv[0]*w[0] + ee[1]*vv[1]*w[1]) * ice_density;

		// Extend layer to proper thickness
		if (V1(i,j) < gcm_thk) {
			double v1old = V1(i,j);
			V1(i,j) = gcm_thk;	// Extend layer to proper thickness
			double fact = V1(i,j) / v1old;
			M1(i,j) *= fact;
			H1(i,j) *= fact;
		}

		// ----- Next layer for GCM
		ierr = integrate_weights_2(vv, gcm_thk, gcm_thk+gcm_thk, w); CHKERRQ(ierr);
		V2(i,j) = vv[0]*w[0] + vv[1]*w[1];
		M2(i,j) = V2(i,j) * ice_density;
		H2(i,j) = (ee[0]*vv[0]*w[0] + ee[1]*vv[1]*w[1]) * ice_density;

		// Extend layer to proper thickness
		if (V2(i,j) == 0) {
			V2(i,j) = V1(i,j);
			M2(i,j) = M1(i,j);
			H2(i,j) = H1(i,j);
		} else if (V2(i,j) < gcm_thk) {
			double v1old = V2(i,j);
			V2(i,j) = gcm_thk;	// Extend layer to proper thickness
			double fact = V2(i,j) / v1old;
			M2(i,j) *= fact;
			H2(i,j) *= fact;
		}

	}}
	ierr = ice_thickness.end_access(); CHKERRQ(ierr);
	ierr = M1.end_access(); CHKERRQ(ierr);
	ierr = M2.end_access(); CHKERRQ(ierr);
	ierr = H1.end_access(); CHKERRQ(ierr);
	ierr = H2.end_access(); CHKERRQ(ierr);
	ierr = V1.end_access(); CHKERRQ(ierr);
	ierr = V2.end_access(); CHKERRQ(ierr);
	ierr = Enth3.end_access(); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode PISMIceModel::grid_setup()
{
	super::grid_setup();

	// super::grid_setup() trashes grid.time->start().  Now set it correctly.
	grid.time->set_start(params.time_start_s);
	grid.time->set(params.time_start_s);

	return 0;
}


PetscErrorCode PISMIceModel::misc_setup()
{
	PetscErrorCode ierr;
	ierr = super::misc_setup(); CHKERRQ(ierr);


	// ------ Initialize MassEnth structures: base, cur, rate
	for (auto &ii : cur.all_vecs) {
		ierr = ii.vec.set(0); CHKERRQ(ierr);
	}
	ierr = compute_enth2(cur.total.enth, cur.total.mass); CHKERRQ(ierr);
	ierr = cur.set_epsilon(grid); CHKERRQ(ierr);

	// base = cur
	auto base_ii(base.all_vecs.begin());
	auto cur_ii(cur.all_vecs.begin());
	for (; base_ii != base.all_vecs.end(); ++base_ii, ++cur_ii) {
		cur_ii->vec.copy_to(base_ii->vec);
	}

	// --------- Set up to write pism_state.nc
	PSConstantGLINT2 *surface = ps_constant_glint2();
	std::vector<IceModelVec *> vecs;
	vecs.push_back(&Enth3);
	vecs.push_back(&ice_thickness);
	vecs.push_back(&ice_surface_temp);
	vecs.push_back(&surface->surface_temp);
	for (auto ii = rate.all_vecs.begin(); ii != rate.all_vecs.end(); ++ii) {
		vecs.push_back(&ii->vec);
	}
	std::string ofname = (params.output_dir / "pism_out_state.nc").string();
	pism_out_state_nc.reset(new VecBundleWriter(&grid, ofname, std::move(vecs)));
	ierr = pism_out_state_nc->init(); CHKERRQ(ierr);


	return 0;
}

/** Sums over columns to compute enthalpy on 2D grid.

NOTE: Unfortunately so far PISM does not keep track of enthalpy in
"partially-filled" cells, so Enth2(i,j) is not valid at locations like
this one. We need to address this, but so far, it seems to me, the
best thing we can do is estimate Enth2(i,j) at partially-filled cells
by computing the average over icy neighbors. I think you can re-use
the idea from IceModel::get_threshold_thickness(...) (iMpartgrid.cc).  */


PetscErrorCode PISMIceModel::compute_enth2(pism::IceModelVec2S &enth2, pism::IceModelVec2S &mass2)
{
	PetscErrorCode ierr;

	//	 getInternalColumn() is allocated already
	double ice_density = config.get("ice_density");
	ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
	ierr = Enth3.begin_access(); CHKERRQ(ierr);
	ierr = enth2.begin_access(); CHKERRQ(ierr);
	ierr = mass2.begin_access(); CHKERRQ(ierr);
	for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
		for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
			enth2(i,j) = 0;
			mass2(i,j) = 0;

			// count all ice, including cells that have so little they
			// are considered "ice-free"
			if (ice_thickness(i,j) > 0) {
				double *Enth;	// do NOT delete this pointer: space returned by
				const int ks = grid.kBelowHeight(ice_thickness(i,j));
				ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
				for (int k=0; k<ks; ++k) {
					double dz = (grid.zlevels[k+1] - grid.zlevels[k]);
					enth2(i,j) += Enth[k] * dz;		// m J / kg
				}

				// Do the last layer a bit differently
				double dz = (ice_thickness(i,j) - grid.zlevels[ks]);
				enth2(i,j) += Enth[ks] * dz;
				enth2(i,j) *= ice_density;		// --> J/m^2
				mass2(i,j) = ice_thickness(i,j) * ice_density;		// --> kg/m^2
			}
		}
	}
	ierr = ice_thickness.end_access(); CHKERRQ(ierr);
	ierr = Enth3.end_access(); CHKERRQ(ierr);
	ierr = enth2.end_access(); CHKERRQ(ierr);
	ierr = mass2.end_access(); CHKERRQ(ierr);
	return 0;
}


/** Merges surface temperature derived from Enth3 into any NaN values
in the vector provided.
@param deltah IN: Input from Glint2 (change in enthalpy of each grid
	cell over the timestep) [W m-2].
@param default_val: The value that deltah(i,j) will have if no value
	is listed for that grid cell
@param timestep_s: Length of the current coupling timestep [s]
@param surface_temp OUT: Resulting surface temperature to use as the Dirichlet B.C.
*/
PetscErrorCode PISMIceModel::construct_surface_temp(
	pism::IceModelVec2S &deltah,			// IN: Input from Glint2
	double default_val,
	double timestep_s,		// Length of this coupling interval [s]
	pism::IceModelVec2S &surface_temp)	// OUT: Temperature @ top of ice sheet (to use for Dirichlet B.C.)

{
	PetscErrorCode ierr;

	printf("BEGIN PISMIceModel::merge_surface_temp default_val=%g\n", default_val);

	PSConstantGLINT2 *surface = ps_constant_glint2();

	double ice_density = config.get("ice_density");

	ierr = Enth3.begin_access(); CHKERRQ(ierr);
	ierr = deltah.begin_access(); CHKERRQ(ierr);
	ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
	ierr = surface_temp.begin_access(); CHKERRQ(ierr);

	// First time around, set effective_surface_temp to top temperature
	for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
	for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
		double &surface_temp_ij(surface_temp(i,j));
		double const &deltah_ij(deltah(i,j));

		double *Enth;
		ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);

		// Enthalpy at top of ice sheet
		const int ks = grid.kBelowHeight(ice_thickness(i,j));
		double spec_enth3 = Enth[ks];		// Specific enthalpy [J kg-1]

		if (deltah_ij != default_val) {
			// Adjust enthalpy @top by deltah
			double toplayer_dz = ice_thickness(i,j) - grid.zlevels[ks];		// [m]

			// [J kg-1] = [J kg-1]
			//     + [J m-2] * [m^2 m-3] * [m^3 kg-1]
			spec_enth3 = spec_enth3
				+ deltah_ij / (toplayer_dz * ice_density) ; //* timestep_s;
		}


		// Convert specific enthalpy value to surface temperature
		double tg3;		// Temperature at top of ice sheet
		const double p = 0.0;		// Pressure at top of ice sheet
		EC->getAbsTemp(spec_enth3, p, surface_temp_ij);
	}}
	ierr = surface_temp.end_access(); CHKERRQ(ierr);
	ierr = ice_thickness.end_access(); CHKERRQ(ierr);
	ierr = deltah.end_access(); CHKERRQ(ierr);
	ierr = Enth3.end_access(); CHKERRQ(ierr);

	printf("END PISMIceModel::merge_surface_temp\n");

	return 0;
}


}}	// namespace glint2::gpism
