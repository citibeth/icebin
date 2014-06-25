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

#define NO_IMPORT_ARRAY
#include "_glint2_module.hpp"

#include "Grid_py.hpp"
#include <structmember.h>	// Also Python-related
#include <numpy/arrayobject.h>
#include <math.h>
#include <glint2/Grid.hpp>
#include "pyutil.hpp"

using namespace glint2;

#define RETURN_INVALID_ARGUMENTS(fname) {\
	fprintf(stderr, fname "(): invalid arguments.\n"); \
	PyErr_SetString(PyExc_ValueError, fname "(): invalid arguments."); \
	return NULL; }


// ========================================================================

void PyGrid::init(std::unique_ptr<glint2::Grid> &&_grid)
{
	grid = std::move(_grid);
//	ncells_full = grid->ncells_full();
}

// ========= class snowdrift.Grid :
static PyObject *Grid_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	PyGrid *self;

	self = new (type->tp_alloc(type, 0)) PyGrid;

//	Py_INCREF(self);
    return (PyObject *)self;
}

/** Read from a file */
static int Grid__init(PyGrid *self, PyObject *args, PyObject *kwds)
{
//printf("Grid__init() called\n");

	// Get arguments
	const char *fname;
	const char *vname;
	if (!PyArg_ParseTuple(args, "ss", &fname, &vname)) {
		// Throw an exception...
		PyErr_SetString(PyExc_ValueError,
			"Grid__init() called without a valid string as argument.");
		return 0;
	}

	// Instantiate pointer
	NcFile nc(fname, NcFile::ReadOnly);
//	std::unique_ptr<glint2::Grid> sptr(read_grid(nc, std::string(vname)).release());
	self->init(read_grid(nc, std::string(vname)));
//fprintf(stderr, "Grid_new() returns %p\n", self->grid);
	nc.close();

	return 0;
}

static void Grid_dealloc(PyGrid *self)
{
	self->~PyGrid();
	self->ob_type->tp_free((PyObject *)self);
}

//static PyMemberDef Grid_members[] = {{NULL}};

// -----------------------------------------------------------

static PyObject * Grid_get_proj_areas(PyGrid *self_py, PyObject *args)
{
	PyObject *ret_py = NULL;
	try {
		// Get Arguments
		char *sproj_py;
		if (!PyArg_ParseTuple(args, "s",
			&sproj_py))
		{
			RETURN_INVALID_ARGUMENTS();
		}

		// Do the call
		auto ret(self_py->grid->get_proj_areas(std::string(sproj_py)));

		// Convert to Python format
		ret_py = giss::vector_to_py(ret);
//printf("Grid_get_proj_areas() finished OK, ret_py = %p\n", ret_py);
		return ret_py;
	} catch(...) {
		printf("Grid_get_proj_areas() encountered an exception...\n");
		if (ret_py) Py_DECREF(ret_py);
		return NULL;
	}
}
// -----------------------------------------------------------
static PyObject * Grid_get_native_areas(PyGrid *self_py, PyObject *args)
{
	PyObject *ret_py = NULL;
	try {
		// Get Arguments
		char *sproj_py;
		if (!PyArg_ParseTuple(args, ""))	// No arguments
		{ return NULL; }

		// Do the call
		auto ret(self_py->grid->get_native_areas());

		// Convert to Python format
		ret_py = giss::vector_to_py(ret);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		return NULL;
	}
}
// -----------------------------------------------------------

static PyMethodDef Grid_methods[] = {
	{"get_proj_areas", (PyCFunction)Grid_get_proj_areas, METH_VARARGS,
		""},
	{"get_native_areas", (PyCFunction)Grid_get_native_areas, METH_VARARGS,
		""},

	{NULL}     /* Sentinel - marks the end of this structure */
};

static PyMemberDef Grid_members[] = {
//	{"ncells_full", T_INT, offsetof(PyGrid, ncells_full)},
	{NULL}
};

static char *Grid_doc =
	"Encapsulates all information about a grid in Glint2.\n"
	"Python wrapper of glint2::Grid.\n"
	"\n"
	"Constructor: Grid(fname, vname)\n"
	"    Read from an existing Grid file.\n"
	"    fname : str\n"
	"        Name of Glint2 configuration file>\n"
	"    vname (OPTIONAL):\n"
	"        Name of variable to read in config file\n";


PyTypeObject GridType = {
   PyObject_HEAD_INIT(NULL)
   0,                         /* ob_size */
   "Grid",               /* tp_name */
   sizeof(PyGrid),     /* tp_basicsize */
   0,                         /* tp_itemsize */
   (destructor)Grid_dealloc, /* tp_dealloc */
   0,                         /* tp_print */
   0,                         /* tp_getattr */
   0,                         /* tp_setattr */
   0,                         /* tp_compare */
   0,                         /* tp_repr */
   0,                         /* tp_as_number */
   0,                         /* tp_as_sequence */
   0,                         /* tp_as_mapping */
   0,                         /* tp_hash */
   0,                         /* tp_call */
   0,                         /* tp_str */
   0,                         /* tp_getattro */
   0,                         /* tp_setattro */
   0,                         /* tp_as_buffer */
   Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags*/
   Grid_doc,                  /* tp_doc */
   0,                         /* tp_traverse */
   0,                         /* tp_clear */
   0,                         /* tp_richcompare */
   0,                         /* tp_weaklistoffset */
   0,                         /* tp_iter */
   0,                         /* tp_iternext */
   Grid_methods,              /* tp_methods */
   Grid_members,              /* tp_members */
//   0,                       /* tp_members */
   0,                         /* tp_getset */
   0,                         /* tp_base */
   0,                         /* tp_dict */
   0,                         /* tp_descr_get */
   0,                         /* tp_descr_set */
   0,                         /* tp_dictoffset */
   (initproc)Grid__init,  /* tp_init */
   0,                         /* tp_alloc */
   (newfunc)&Grid_new    /* tp_new */
//   (freefunc)Grid_free	/* tp_free */
};
