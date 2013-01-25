#include <vector>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <cmath>
#include <giss/pyutil.hpp>

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

// All the classes (PyTypeObject) for the entire module
static PyTypeObject *glint2_types[] = {
	NULL
};

// ===========================================================================

extern "C"
void initglint2(void)
{
	giss::init_module("_glint2",
		"Interface to C++ GLINT2 Library",
		glint2_function_sets, glint2_types);
}
