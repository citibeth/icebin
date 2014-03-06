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
namespace pism {


/** This is the GLINT2 customized version of PISM's ::IceModel class.

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

class PISMIceModel : public ::IceModel
{
public:

	// see iceModel.cc for implementation of constructor and destructor:
	PISMIceModel(IceGrid &g, PISMConfig &config, PISMConfig &overrides);
	virtual ~PISMIceModel(); // must be virtual merely because some members are virtual

	virtual PetscErrorCode allocate_subglacial_hydrology();
	virtual PetscErrorCode allocate_couplers();

	/** @return Our instance of PSConstantGLINT2 */
	PSConstantGLINT2 *ps_constant_glint2()
		{ return dynamic_cast<PSConstantGLINT2 *>(surface); }
	NullTransportHydrology *null_hydrology()
		{ return dynamic_cast<NullTransportHydrology *>(::IceModel::subglacial_hydrology); }
};

}}
