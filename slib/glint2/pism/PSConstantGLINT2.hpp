// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef _PSCONSTANTGLINT2_H_
#define _PSCONSTANTGLINT2_H_

#include <PISMSurface.hh>
#include <iceModelVec.hh>
#include <PISMAtmosphere.hh>


namespace glint2 {
namespace gpism {

//! \brief A class implementing a constant-in-time surface model for the surface mass balance.
//!
//! Reads data from a PISM input file.
//!
//! Ice surface temperature is parameterized as in PISM-GLINT2, using a latitude
//! and surface elevation-dependent formula.

class PSConstantGLINT2 : public pism::SurfaceModel {
	bool _initialized;

public:
	/** Experiment with b.c. in terms of conductive heat flow. */
//	virtual BCType get_conduction_bc_type() { return NEUMANN; }
	virtual BCType get_conduction_bc_type() { return DIRICHLET; }

	/** @param conf Not Used (Looked up all the constructors, it just
	sets this->config, whic his not used */
	PSConstantGLINT2(pism::IceGrid &g, const pism::Config &conf);

	virtual PetscErrorCode init(pism::Vars &vars);

	/** Just deletes input, no atmosphere is needed for this surface model. */
	virtual void attach_atmosphere_model(pism::AtmosphereModel *input);

	virtual void get_diagnostics(std::map<std::string, pism::Diagnostic*> &dict,
															 std::map<std::string, pism::TSDiagnostic*> &ts_dict);
	virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);
	virtual PetscErrorCode ice_surface_mass_flux(pism::IceModelVec2S &result);
//	virtual PetscErrorCode ice_surface_heat_flux(pism::IceModelVec2S &result);
	virtual PetscErrorCode ice_surface_temperature(pism::IceModelVec2S &result);
	virtual PetscErrorCode define_variables(std::set<std::string> vars, const pism::PIO &nc, pism::IO_Type nctype);
	virtual PetscErrorCode write_variables(std::set<std::string> vars, const pism::PIO &nc);
	virtual void add_vars_to_output(std::string keyword, std::set<std::string> &result);
protected:
	std::string input_file;
public:
	// Inputs from Glint2
	pism::IceModelVec2S glint2_smb_mass;
	pism::IceModelVec2S glint2_surface_temp;
	pism::IceModelVec2S glint2_heat_flux;

	// IMPLIED: liquid fraction of 0 (see our superclass)
	pism::IceModelVec2S _ice_surface_hflux;
private:
	PetscErrorCode allocate_PSConstantGLINT2();
};

}}		// namespace

#endif /* _PSCONSTANTGLINT2_H_ */
