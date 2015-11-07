#pragma once

// --------------------------------
// PISM Includes... want to be included first
#include <petsc.h>
#include "IceGrid.hh"
#include "PISMNCFile.hh"
#include "iceModel.hh"
// --------------------------------

#include <vector>
#include <string>

namespace glint2 {
namespace gpism {


/** Sets up to easily write out a bundle of PISM variables to a file. */
class VecBundleWriter {
	pism::IceGrid *grid;
	std::string fname;			// Name of the file to write
	std::vector<pism::IceModelVec *> vecs;	// The vectors we will write

public:

	VecBundleWriter(
		pism::IceGrid *_grid,
		std::string const &_fname,
		std::vector<pism::IceModelVec *> &&_vecs);

	PetscErrorCode init();

	/** Dump the value of the Vectors at curent PISM simulation time. */
	PetscErrorCode write();
};


}}
