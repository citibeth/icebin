/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "_glint2_module.hpp"

#include <vector>
#include <cmath>
#include "pyutil.hpp"

// ============================================================
// Type definitions for the Python classes this module defines

// Function defintions, per sub-module
extern PyMethodDef matrix_ops_functions[];

// ===========================================================================

// All the functions for the entire module
static PyMethodDef *glint2_function_sets[] = {
	matrix_ops_functions,
	NULL
};

extern PyTypeObject NcFileType;
extern PyTypeObject GridType;
extern PyTypeObject MatrixMakerType;

// All the classes (PyTypeObject) for the entire module
static PyTypeObject *glint2_types[] = {
	&NcFileType,
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
