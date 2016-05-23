#pragma once

#include <petsc.h>
#include <mpi.h>

namespace glint2 {
namespace gpism {


/* This explicit scoping forces destructors to be called before PetscFinalize() */
class PetscContext {
protected:
	MPI_Comm _comm;

public:
	PetscContext(MPI_Comm comm, int &argc, char**&argv) : _comm(comm) {
		if (allocate(argc, argv) != 0) {
			PetscPrintf(_comm, "PetscContext::PetscContext(...): allocate() failed\n");
// Don't know how we should tie together error reporting between Petsc/PISM, ModelE and GLINT2
//			PISMEnd();
		}
	}

	~PetscContext() {
		if (deallocate() != 0) {
			PetscPrintf(_comm, "PetscContext::~PetscContext(...): deallocate() failed\n");
//			PISMEnd();
		}
	}

	int allocate(int &argc, char **&argv) {
		static char help[] =
			"GLINT2 PISM module";
		PetscErrorCode ierr;
printf("Calling PetscInitialize\n");
#if 0
        /* See note in IceModel_PISM.cpp for newer versions of PETSc... */
        ierr = PetscOptionsSetValue("-no_signal_handler","true");
#endif
		ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);
printf("Done calling PetscInitialize\n");
		return 0;
	}

	int deallocate() {
		PetscErrorCode ierr;
		ierr = PetscFinalize(); CHKERRQ(ierr);
		return 0;
	}
};

}}
