// Copyright (C) 2008-2014 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.	See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA	02110-1301	USA

#include "PSConstantGLINT2.hpp"
#include <PIO.hh>
#include <PISMVars.hh>
#include <IceGrid.hh>

using namespace pism;

namespace glint2 {
namespace gpism {

PSConstantGLINT2::PSConstantGLINT2(pism::IceGrid &g, const pism::Config &conf)
	: pism::SurfaceModel(g, conf), _initialized(false)
{
	PetscErrorCode ierr = allocate_PSConstantGLINT2(); CHKERRCONTINUE(ierr);
	if (ierr != 0) {
		PISMEnd();
	}
}

void PSConstantGLINT2::attach_atmosphere_model(pism::AtmosphereModel *input)
{
	delete input;
}

PetscErrorCode PSConstantGLINT2::allocate_PSConstantGLINT2()
{
printf("BEGIN PSConstantGLINT2::allocate_PSConstantGLINT2()\n");
	PetscErrorCode ierr;

printf("PSConstantGLINT2::allocate(): grid=%p, Mx My = %d %d\n", &grid, grid.Mx, grid.My);

	ierr = glint2_deltah.create(grid, "glint2_deltah", WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = glint2_deltah.set_attrs("climate_state",
		"enthalpy of constant-in-time ice-equivalent surface mass balance (accumulation/ablation) rate",
		"J m-2", ""); CHKERRQ(ierr);
	ierr = glint2_deltah.set_glaciological_units("kg m-2 year-1"); CHKERRQ(ierr);
//	glint2_deltah.write_in_glaciological_units = true;


	ierr = glint2_massxfer_rate.create(grid, "glint2_massxfer_rate", WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = glint2_massxfer_rate.set_attrs("climate_state",
		"enthalpy of constant-in-time ice-equivalent surface mass balance (accumulation/ablation) rate",
		"kg m-2 s-1", ""); CHKERRQ(ierr);


	ierr = glint2_enthxfer_rate.create(grid, "glint2_enthxfer_rate", WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = glint2_enthxfer_rate.set_attrs("climate_state",
		"constant-in-time heat flux through top surface",
		"W m-2", ""); CHKERRQ(ierr);

	// This variable is computed from the inputs above.
	ierr = surface_temp.create(grid, "surface_temp", WITHOUT_GHOSTS); CHKERRQ(ierr);
	ierr = surface_temp.set_attrs("climate_state",
		"Temperature to use for Dirichlet B.C. at surface",
		"K", ""); CHKERRQ(ierr);

printf("END PSConstantGLINT2::allocate_PSConstantGLINT2()\n");
	return 0;
}

PetscErrorCode PSConstantGLINT2::init(pism::Vars &vars)
{
	// This is called (via pism::IceModel::init_couplers()) from both
	// pism::IceModel::misc_setup() and pism::IceModel::model_state_setup()
	if (_initialized) return 0;

printf("BEGIN PSConstantGLINT2::init(this=%p)\n", this);
	PetscErrorCode ierr;
	bool do_regrid = false;
	int start = -1;

	m_t = m_dt = GSL_NAN;	// every re-init restarts the clock

	ierr = verbPrintf(2, grid.com,
		 "* Initializing the PSConstantGLINT2 surface model. Serves as storage for climate fields.\n"
		 "	Any choice of atmosphere coupler (option '-atmosphere') is ignored.\n"); CHKERRQ(ierr);

	// find PISM input file to read data from:
	// (do_regrid will be true if there's a -boot_file
	// command line option.  We are not using it that way with PISM)
	// see pism/src/base/util/PISMComponent.cc
	ierr = find_pism_input(input_file, do_regrid, start); CHKERRQ(ierr);

	// It doesn't matter what we set this to, it will be re-set later.
	ierr = glint2_deltah.set(0.0); CHKERRQ(ierr);
	ierr = glint2_massxfer_rate.set(0.0); CHKERRQ(ierr);
	ierr = glint2_enthxfer_rate.set(0.0); CHKERRQ(ierr);
	ierr = surface_temp.set(0.0); CHKERRQ(ierr);

	// parameterizing the ice surface temperature 'ice_surface_temp'
	ierr = verbPrintf(2, grid.com,
				"		parameterizing the ice surface temperature 'ice_surface_temp' ... \n"); CHKERRQ(ierr);

printf("END PSConstantGLINT2::init()\n");
	_initialized = true;
	return 0;
}

PetscErrorCode PSConstantGLINT2::update(PetscReal my_t, PetscReal my_dt)
{
//	PetscErrorCode ierr;

	if ((fabs(my_t - m_t) < 1e-12) &&
			(fabs(my_dt - m_dt) < 1e-12))
		return 0;

	m_t	= my_t;
	m_dt = my_dt;

	return 0;
}

void PSConstantGLINT2::get_diagnostics(std::map<std::string, pism::Diagnostic*> &/*dict*/,
	std::map<std::string, pism::TSDiagnostic*> &/*ts_dict*/)
{
	// empty (does not have an atmosphere model)
}

  // Returns [kg m-2 s-1]
PetscErrorCode PSConstantGLINT2::ice_surface_mass_flux(IceModelVec2S &result) {
	PetscErrorCode ierr;

	ierr = glint2_massxfer_rate.copy_to(result); CHKERRQ(ierr);
	return 0;
}

PetscErrorCode PSConstantGLINT2::ice_surface_temperature(IceModelVec2S &result) {
	PetscErrorCode ierr;

	ierr = surface_temp.copy_to(result); CHKERRQ(ierr);
	return 0;
}

// PetscErrorCode PSConstantGLINT2::ice_surface_heat_flux(IceModelVec2S &result) {
// 	PetscErrorCode ierr;
// 
// 	ierr = glint2_enthxfer_rate.copy_to(result); CHKERRQ(ierr);
// 	return 0;
// }

void PSConstantGLINT2::add_vars_to_output(std::string /*keyword*/, std::set<std::string> &result) {
	result.insert("glint2_deltah");
	result.insert("glint2_massxfer_rate");
	result.insert("glint2_enthxfer_rate");
	result.insert("surface_temp");
	// does not call atmosphere->add_vars_to_output().
}

PetscErrorCode PSConstantGLINT2::define_variables(std::set<std::string> vars, const PIO &nc, pism::IO_Type nctype) {
	PetscErrorCode ierr;

	ierr = pism::SurfaceModel::define_variables(vars, nc, nctype); CHKERRQ(ierr);

	if (set_contains(vars, "glint2_enthxfer_rate")) {
		ierr = glint2_enthxfer_rate.define(nc, nctype); CHKERRQ(ierr);
	}

	if (set_contains(vars, "glint2_deltah")) {
		ierr = glint2_deltah.define(nc, nctype); CHKERRQ(ierr);
	}

	if (set_contains(vars, "glint2_massxfer_rate")) {
		ierr = glint2_massxfer_rate.define(nc, nctype); CHKERRQ(ierr);
	}

	if (set_contains(vars, "surface_temp")) {
		ierr = surface_temp.define(nc, nctype); CHKERRQ(ierr);
	}
	return 0;
}

PetscErrorCode PSConstantGLINT2::write_variables(std::set<std::string> vars, const PIO &nc) {
	PetscErrorCode ierr;

	if (set_contains(vars, "glint2_enthxfer_rate")) {
		ierr = glint2_enthxfer_rate.write(nc); CHKERRQ(ierr);
	}

	if (set_contains(vars, "glint2_deltah")) {
		ierr = glint2_deltah.write(nc); CHKERRQ(ierr);
	}

	if (set_contains(vars, "glint2_massxfer_rate")) {
		ierr = glint2_massxfer_rate.write(nc); CHKERRQ(ierr);
	}

	if (set_contains(vars, "surface_temp")) {
		ierr = surface_temp.write(nc); CHKERRQ(ierr);
	}

	return 0;
}

}		// namespace glint2::gpism
}		// namespace glint2
