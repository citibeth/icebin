#include <glint2/pism/PISMIceModel.hpp>

namespace glint2 {
namespace pism {

PISMIceModel::PISMIceModel(IceGrid &g, PISMConfig &config, PISMConfig &overrides) :
	::IceModel(g, config, overrides)
{}

PISMIceModel::~PISMIceModel() {} // must be virtual merely because some members are virtual


PetscErrorCode PISMIceModel::allocate_subglacial_hydrology()
{
	printf("BEGIN PISMIceModel::allocate_subglacial_hydrology()\n");
//	printf("subglacial_hydrology = %p %p\n", subglacial_hydrology, ::IceModel::subglacial_hydrology);
	if (::IceModel::subglacial_hydrology != NULL) return 0; // indicates it has already been allocated
    subglacial_hydrology = new NullTransportHydrology(grid, config);
	printf("END PISMIceModel::allocate_subglacial_hydrology()\n");
	return 0;
}



PetscErrorCode PISMIceModel::allocate_couplers()
{
	PetscErrorCode ierr;
	// Initialize boundary models:
	PAFactory pa(grid, config);
	PSFactory ps(grid, config);
	POFactory po(grid, config);
	PISMAtmosphereModel *atmosphere;

	ierr = PetscOptionsBegin(grid.com, "", "Options choosing PISM boundary models", ""); CHKERRQ(ierr);

#if 1
	// GLINT2-modified version
	if (surface == NULL) {
		surface = new PSConstantGLINT2(grid, config);
		external_surface_model = false;

		pa.create(atmosphere);
		surface->attach_atmosphere_model(atmosphere);
	}
#else
	// Original Version
	if (surface == NULL) {
		ps.create(surface);
		external_surface_model = false;

		pa.create(atmosphere);
		surface->attach_atmosphere_model(atmosphere);
	}
#endif

	if (ocean == NULL) {
		po.create(ocean);
		external_ocean_model = false;
	}
	ierr = PetscOptionsEnd(); CHKERRQ(ierr);

	return 0;
}

}}
