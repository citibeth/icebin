#include "VecBundleWriter.hpp"
#include "PISMTime.hh"

using namespace pism;

namespace glint2 {
namespace gpism {

VecBundleWriter::VecBundleWriter(
	pism::IceGrid *_grid,
	std::string const &_fname,
	std::vector<pism::IceModelVec *> &&_vecs)
: grid(_grid), fname(_fname), vecs(std::move(_vecs))
{}

PetscErrorCode VecBundleWriter::init()
{
	PetscErrorCode ierr;
	pism::PIO nc(*grid, grid->config.get_string("output_format"));

	ierr = nc.open(fname, PISM_READWRITE_MOVE); CHKERRQ(ierr);
    ierr = nc.def_time(grid->config.get_string("time_dimension_name"),
		grid->time->calendar(),
		grid->time->CF_units_string()); CHKERRQ(ierr);

	for (pism::IceModelVec *vec : vecs) {
		ierr = vec->define(nc, PISM_DOUBLE); CHKERRQ(ierr);
	}

	// --------- Close and return
	nc.close();

	return 0;
}

/** Dump the value of the Vectors at curent PISM simulation time. */
PetscErrorCode VecBundleWriter::write()
{
	PetscErrorCode ierr;
	double t1;
	pism::PIO nc(*grid, grid->config.get_string("output_format"));

	ierr = nc.open(fname.c_str(), PISM_READWRITE); CHKERRQ(ierr); // append to file
	ierr = nc.append_time(grid->config.get_string("time_dimension_name"), t1); CHKERRQ(ierr);

	for (pism::IceModelVec *vec : vecs) {
		ierr = vec->write(nc, PISM_DOUBLE); CHKERRQ(ierr);
	}
	nc.close();

	return 0;
}

}}
