#include "_glint2_module.hpp"

#include <vector>
#include <cmath>
#include "pyutil.hpp"

// ============================================================
// Type definitions for the Python classes this module defines
// extern PyTypeObject GridType

// Function defintions, per sub-module
extern PyMethodDef matrix_ops_functions[];

// ===========================================================================

// All the functions for the entire module
static PyMethodDef *glint2_function_sets[] = {
	matrix_ops_functions,
	NULL
};

extern PyTypeObject GridType;
extern PyTypeObject MatrixMakerType;

// All the classes (PyTypeObject) for the entire module
static PyTypeObject *glint2_types[] = {
	&GridType,
	&MatrixMakerType,
	NULL
};

// ===========================================================================

extern "C"
void init_glint2(void)
{
	giss::init_module("_glint2",
		"Interface to C++ GLINT2 Library",
		glint2_function_sets, glint2_types);


	/* See http://dsnra.jpl.nasa.gov/software/Python/numpydoc/numpy-13.html

	In addition to including arrayobject.h , the extension must call
	import_array() in its initialization function, after the call to
	Py_InitModule() . This call makes sure that the module which
	implements the array type has been imported, and initializes a pointer
	array through which the NumPy functions are called. If you forget this
	call, your extension module will crash on the first call to a NumPy
	function. */

	import_array();
}
