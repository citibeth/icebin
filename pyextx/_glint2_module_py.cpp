/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
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

#include <cstdio>
#include <vector>
#include <cmath>
#include "pyutil.hpp"
#include <giss/exit.hpp>
#include <everytrace.h>

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
// Ported to Python3
// https://docs.python.org/3/howto/cporting.html


struct module_state {
    PyObject *error;
};
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

static int glint2Traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int glint2Clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static PyModuleDef glint2ModuleDef = {
    PyModuleDef_HEAD_INIT,
    "_glint2",
    "Interface to C++ GLINT2 Library",
    sizeof(struct module_state),
    NULL,
    NULL,
    glint2Traverse,
    glint2Clear,
    NULL
};

// ===========================================================================

/** Exit subroutine to use when we're running Python */
static void python_exit(int i)
{
    everytrace_dump();
    fprintf(stderr, "Returning to Python interpreter.\n");
    throw std::exception();     // This will get caught, then go back to Python
}


extern "C"
PyObject *PyInit__glint2(void)
{
//  libglint2_ncerror_segfault();

    giss::exit = &python_exit;


    PyObject *mod = giss::init_module(glint2ModuleDef,
        glint2_function_sets, glint2_types);

    /* See http://dsnra.jpl.nasa.gov/software/Python/numpydoc/numpy-13.html

    In addition to including arrayobject.h , the extension must call
    import_array() in its initialization function, after the call to
    Py_InitModule() . This call makes sure that the module that
    implements the array type has been imported, and initializes a pointer
    array through which the NumPy functions are called. If you forget this
    call, your extension module will crash on the first call to a NumPy
    function. */
    import_array();

    return (PyObject *)mod;
}

