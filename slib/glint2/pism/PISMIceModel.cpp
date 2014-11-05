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

#if 1
//Temporarily comment out.  For now, this creates problems in ncView.
//strain_heating is inf at the coastlines.
// See: https://github.com/pism/pism/issues/292
	// ------------ Volumetric Strain Heating
	// strain_heating_sum += my_dt * sum_columns(strainheating3p)
	IceModelVec3 *strain_heating3p;
	stress_balance->get_volumetric_strain_heating(strain_heating3p);
	// cur.strain_heating = cur.strain_heating * 1.0 + my_dt * sum_columns(strain_heating3p)
	ierr = strain_heating3p->sumColumns(cur.strain_heating, 1.0, my_dt); CHKERRQ(ierr);
#endif

	printf("END PISMIceModel::energyStep(time=%f)\n", t_TempAge);
	return 0;
}

PetscErrorCode PISMIceModel::massContExplicitStep() {
	PetscErrorCode ierr;

	printf("BEGIN PISMIceModel::MassContExplicitStep()\n");

	_ice_density = config.get("ice_density");
	_meter_per_s_to_kg_per_m2 = dt * _ice_density;
printf("_meter_per_s_to_kg_per_m2: %g * %g = %g\n", dt, _ice_density, _meter_per_s_to_kg_per_m2);


	// =========== The Mass Continuity Step Itself
	// This will call through to accumulateFluxes_massContExplicitStep()
	// in the inner loop
	ierr = Enth3.begin_access(); CHKERRQ(ierr);
	ierr = cur.pism_smb.begin_access(); CHKERRQ(ierr);
	ierr = cur.melt_grounded.begin_access(); CHKERRQ(ierr);
	ierr = cur.internal_advection.begin_access(); CHKERRQ(ierr);
	ierr = cur.href_to_h.begin_access(); CHKERRQ(ierr);
	ierr = cur.nonneg_rule.begin_access(); CHKERRQ(ierr);

	ierr = super::massContExplicitStep(); CHKERRQ(ierr);

	ierr = cur.nonneg_rule.end_access(); CHKERRQ(ierr);
	ierr = cur.href_to_h.end_access(); CHKERRQ(ierr);
	ierr = cur.internal_advection.end_access(); CHKERRQ(ierr);
	ierr = cur.melt_grounded.end_access(); CHKERRQ(ierr);
	ierr = cur.pism_smb.end_access(); CHKERRQ(ierr);
	ierr = Enth3.end_access(); CHKERRQ(ierr);



	// =========== AFTER the Mass Continuity Step

	// ----------- SMB: mass and enthalpy
	PSConstantGLINT2 *surface = ps_constant_glint2();
	double p_air = EC->getPressureFromDepth(0.0);
	ierr = surface->climatic_mass_balance.begin_access(); CHKERRQ(ierr);
	ierr = surface->ice_surface_temp.begin_access(); CHKERRQ(ierr);
	ierr = cur.surface_mass_balance.begin_access(); CHKERRQ(ierr);
	for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
	for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
		double mass = surface->climatic_mass_balance(i,j);		// Our input is in [kg m-2 s-1]
		double specific_enth;
		double T = surface->ice_surface_temp(i,j);
		ierr = EC->getEnthPermissive(T,
			0.0, p_air, specific_enth);	CHKERRQ(ierr); // [J kg-1]

		cur.surface_mass_balance.mass(i,j) += dt * mass;					// [kg m-2]
		cur.surface_mass_balance.enth(i,j) += dt * mass * specific_enth;	// [J m-2]
	}}
	ierr = cur.surface_mass_balance.end_access(); CHKERRQ(ierr);
	ierr = surface->ice_surface_temp.end_access(); CHKERRQ(ierr);
	ierr = surface->climatic_mass_balance.end_access(); CHKERRQ(ierr);

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

	// -------------- Get the easy veriables out of the way...
	cur.pism_smb(i,j) += surface_mass_balance * _meter_per_s_to_kg_per_m2;
	cur.nonneg_rule(i,j) -= nonneg_rule_flux * _ice_density;
	cur.href_to_h(i,j) += Href_to_H_flux * _ice_density;


	// -------------- Melting
	double p_basal = EC->getPressureFromDepth(ice_thickness(i,j));
	double T = EC->getMeltingTemp(p_basal);
	double specific_enth;
	ierr = EC->getEnthPermissive(T, 1.0, p_basal, specific_enth); CHKERRQ(ierr);
	double mass;

	// ------- Melting at base of ice sheet
	mass = -meltrate_grounded * _meter_per_s_to_kg_per_m2;
	cur.melt_grounded.mass(i,j) += mass;
	cur.melt_grounded.enth(i,j) += mass * specific_enth;

	// ------- Melting under ice shelf
	mass = -meltrate_floating * _meter_per_s_to_kg_per_m2;
	cur.melt_floating.mass(i,j) += mass;
	cur.melt_floating.enth(i,j) += mass * specific_enth;


	// -------------- internal_advection
	const int ks = grid.kBelowHeight(ice_thickness(i,j));
	double *Enth;
	ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
	specific_enth = Enth[ks];		// Approximate, we will use the enthalpy of the top layer...

	mass = -(divQ_SIA + divQ_SSA) * _meter_per_s_to_kg_per_m2;

	cur.internal_advection.mass(i,j) += mass;
	cur.internal_advection.enth(i,j) += mass * specific_enth;


//	printf("END PISMIceModel::accumulateFluxes_MassContExplicitStep()\n");
	return 0;
}



PetscErrorCode PISMIceModel::prepare_nc(std::string const &fname, std::unique_ptr<PIO> &nc)
{
	PetscErrorCode ierr;

	nc.reset(new PIO(grid, grid.config.get_string("output_format")));

	ierr = nc->open(fname, PISM_READWRITE_MOVE); CHKERRQ(ierr);
    ierr = nc->def_time(config.get_string("time_dimension_name"),
                       grid.time->calendar(),
                       grid.time->CF_units_string()); CHKERRQ(ierr);
// These are in iMtimseries, but not listed as required in iceModelVec.hh
//    ierr = nc->put_att_text(config.get_string("time_dimension_name"),
//                           "bounds", "time_bounds"); CHKERRQ(ierr);
//    ierr = write_metadata(nc, true, false); CHKERRQ(ierr);
//	ierr = nc->close(): CHKERRQ(ierr);

	return 0;
}


/** @param t0 Time of last time we coupled. */
PetscErrorCode PISMIceModel::set_rate(double dt)
{
	PetscErrorCode ierr;

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

			// base = cur: For ALL vectors
			vbase(i,j) = vcur(i,j);
		}}
		ierr = vrate.end_access(); CHKERRQ(ierr);
		ierr = vcur.end_access(); CHKERRQ(ierr);
		ierr = vbase.end_access(); CHKERRQ(ierr);
	}

	return 0;
}


/** @param t0 Time of last time we coupled. */
PetscErrorCode PISMIceModel::prepare_outputs(double t0)
{
	PetscErrorCode ierr;

	printf("BEGIN PISMIceModel::prepare_outputs()\n");

	// ------ Difference between now and the last time we were called
	double t1 = enthalpy_t();	// Current time of the enthalpy portion of ice model.
	set_rate(t1 - t0);

	// ====================== Compute "custom" variables in output
	// --------- basal_runoff = melt_grounded + melt_floating
	ierr = rate.melt_grounded.begin_access(); CHKERRQ(ierr);
	ierr = rate.melt_floating.begin_access(); CHKERRQ(ierr);
	ierr = basal_runoff.begin_access(); CHKERRQ(ierr);
	for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
	for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
		basal_runoff.mass(i,j) = rate.melt_grounded.mass(i,j) + rate.melt_floating.mass(i,j);
		basal_runoff.enth(i,j) = rate.melt_grounded.enth(i,j) + rate.melt_floating.enth(i,j);
	}}
	ierr = basal_runoff.end_access(); CHKERRQ(ierr);
	ierr = rate.melt_floating.end_access(); CHKERRQ(ierr);
	ierr = rate.melt_grounded.end_access(); CHKERRQ(ierr);

	// --------- ice_surface_enth from Enth3
	ierr = Enth3.begin_access(); CHKERRQ(ierr);
	ierr = ice_surface_enth.begin_access(); CHKERRQ(ierr);
	ierr = ice_surface_enth_depth.begin_access(); CHKERRQ(ierr);
	ierr = ice_thickness.begin_access(); CHKERRQ(ierr);
	for (int i = grid.xs; i < grid.xs + grid.xm; ++i) {
	for (int j = grid.ys; j < grid.ys + grid.ym; ++j) {
		double *Enth;
		ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);

		const int ks = grid.kBelowHeight(ice_thickness(i,j));
		ice_surface_enth(i,j) = Enth[ks];
		// Depth at which enthalpy is as above (could be zero)
		ice_surface_enth_depth(i,j) = ice_thickness(i,j) - grid.zlevels[ks];
	}}
	ierr = ice_thickness.end_access(); CHKERRQ(ierr);
	ierr = ice_surface_enth_depth.end_access(); CHKERRQ(ierr);
	ierr = ice_surface_enth.end_access(); CHKERRQ(ierr);
	ierr = Enth3.end_access(); CHKERRQ(ierr);

	// ====================== Write to the post_energy.nc file (OPTIONAL)

	// ------ Write it out
	PIO nc(grid, grid.config.get_string("output_format"));
	nc.open((params.output_dir / "post_energy.nc").c_str(), PISM_READWRITE);	// append to file
	nc.append_time(config.get_string("time_dimension_name"), t1);
	Enth3.write(nc, PISM_DOUBLE);
	ice_thickness.write(nc, PISM_DOUBLE);
	ice_surface_temp.write(nc, PISM_DOUBLE);
	PSConstantGLINT2 *surface = ps_constant_glint2();
	surface->ice_surface_temp.write(nc, PISM_DOUBLE);
	for (auto ii = rate.all_vecs.begin(); ii != rate.all_vecs.end(); ++ii) {
		ierr = ii->vec.write(nc, PISM_DOUBLE); CHKERRQ(ierr);
	}
	nc.close();

	printf("END PISMIceModel::prepare_outputs()\n");
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


PetscErrorCode PISMIceModel::allocate_internal_objects()
{
	PetscErrorCode ierr;

	ierr = super::allocate_internal_objects(); CHKERRQ(ierr);


	ierr = base.create(grid, "", pism::WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = cur.create(grid, "", pism::WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = rate.create(grid, "", pism::WITHOUT_GHOSTS); CHKERRQ(ierr);

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
	auto bii(base.all_vecs.begin());
	auto cii(cur.all_vecs.begin());
	for (; bii != base.all_vecs.end(); ++bii, ++cii) {
		cii->vec.copy_to(bii->vec);
	}

	// ---------- Create the netCDF output file
	std::unique_ptr<PIO> nc;
	std::string ofname = (params.output_dir / "post_energy.nc").string();
	ierr = prepare_nc(ofname, nc); CHKERRQ(ierr);

	// -------- Define MethEnth structres in netCDF file
	Enth3.define(*nc, PISM_DOUBLE);
	ice_thickness.define(*nc, PISM_DOUBLE);
	ice_surface_temp.define(*nc, PISM_DOUBLE);
	PSConstantGLINT2 *surface = ps_constant_glint2();
	surface->ice_surface_temp.define(*nc, PISM_DOUBLE);
	for (auto ii = rate.all_vecs.begin(); ii != rate.all_vecs.end(); ++ii) {
		ierr = ii->vec.define(*nc, PISM_DOUBLE); CHKERRQ(ierr);
	}

	// --------- Close and return
	nc->close();
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
//	ierr = ice_surface_temp.begin_access(); CHKERRQ(ierr);
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

}}	// namespace glint2::gpism



// This subroutine might be useful in the future.
// It is like an optimized version of calling compute_enth2() twice.
// /** Computes enth(vH1 configuration) - enth(vH0 configuration).
// vH1 is new height of ice sheet, vH0 is old height.
// NOTE: Caller MUST have already called Enth3.begin_access()!
// @param i The gridcell for which enthalpy difference is to be computed.
// @param j
// @param vH0 Old ice thickness.
// @param vH1 New ice thickness.
// @param dMass OUT: Difference in mass [kg m-2]
// @param dEnth OUT: Difference in enthalpy [J m-2] */
// PetscErrorCode PISMIceModel::enth_diff(int i, int j, double vH0, double vH1,
// double &dMass, double &dEnth)
// {
//   PetscErrorCode ierr;
//   
//   double sign = 1.0;
//   if (vH1 < vH0) {
//     std::swap(vH0, vH1);
//     sign = -1.0;
//   }
// 
//   // This is inefficient...
//   const int ks0 = grid.kBelowHeight(vH0);
//   const int ks1 = grid.kBelowHeight(vH1);
// 
//   double *Enth;  // J/kg do NOT delete this pointer: space returned by
//   ierr = Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
// 
//   if (ks0 == ks1) {
//     return Enth(ks0) * (vH1 - vH0);
//   } else {
//     double sum_enth = 0;
// 
//     // First level is partial
//     zhi = grid.zlevels[ks0+1];
//     sum_enth += Enth(ks0) * (zhi - vH0);
// 
//     // Middle full levels
//     for (ks=ks0+1; ks<ks1; ++ks) {
//       zlo = grid.zlevels[ks0];
//       zhi = grid.zlevels[ks0+1];
//       sum_enth += Enth(ks) * (zhi - zlo);
//     }
// 
//     // Last level is partial
//     zlo = grid.zlevels[ks1];
//     sum_enth += Enth(ks0) * (vH1 - zlo);
// 
//   }
// 
//   dEnth = sum_enth * sign;
//   dMass = (vH1 - vH0) * sign * ice_density;
// 
//   return 0;
// }


