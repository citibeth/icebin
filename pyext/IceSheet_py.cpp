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

// Currently not compiled as part of the extension (Oct 2015).

#define NO_IMPORT_ARRAY
#include "_glint2_module.hpp"

#include "IceSheet_py.hpp"
#include <structmember.h>	// Also Python-related
#include <numpy/arrayobject.h>
#include <math.h>
#include <glint2/IceSheet.hpp>
#include "pyutil.hpp"

using namespace glint2;

#define RETURN_INVALID_ARGUMENTS(fname) {\
	fprintf(stderr, fname "(): invalid arguments.\n"); \
	PyErr_SetString(PyExc_ValueError, fname "(): invalid arguments."); \
	return NULL; }


// ========================================================================

void PyIceSheet::init(std::unique_ptr<glint2::IceSheet> &&_sheet)
{
	sheet = std::move(_sheet);
}

// ========= class snowdrift.IceSheet :
static PyObject *IceSheet_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	PyIceSheet *self;

	self = new (type->tp_alloc(type, 0)) PyIceSheet;

//	Py_INCREF(self);
    return (PyObject *)self;
}

/** Read from a file */
static int IceSheet__init(PyIceSheet *self, PyObject *args, PyObject *kwds)
{
//printf("IceSheet__init() called\n");

	// Get arguments
	const char *grid2_fname_py;
	const char *exgrid_fname_py;
	PyObject *elev2_py = NULL;
	PyObject *mask2_py = NULL;

	if (!PyArg_ParseTuple(args, "ssOO",
	   &grid2_fname_py, &exgrid_fname_py, &elev2_py, &mask2_py))
	{
		// Throw an exception...
		PyErr_SetString(PyExc_ValueError,
			"IceSheet__init() called without valid arguments.");
		return 0;
	}


	blitz::Array<int,1> mask2;
	blitz::Array<int,1> *mask2p = NULL;
	if (mask2_py != Py_None) {
		mask2.reference(giss::py_to_blitz<int,1>(mask2_py, "mask2", {overlap.ncol}));
		mask2p = &mask2;
	}

	// Instantiate C++ Ice Sheet
	std::unique_ptr<glint2::IceSheet> sheet(new IceSheet);

	// Fill it in...
	sheet->grid2 = read_grid(NcFile(grid2_fname, NcFile::ReadOnly), "grid");
	sheet->exgrid = read_grid(NcFile(exgrid_fname, NcFile::ReadOnly), "grid");

	int n2 = sheet->grid2->ndata();
	if (mask2_py != Py_None) {
		sheet->mask2.reset(
			new blitz::Array<int,1>(
				giss::py_to_blitz<int,1>(mask2_py, "mask2", {n2}).copy()
			));
	}

	sheet->elev2 = giss::py_to_blitz<double,1>(
		elev2_py, "elev2", {n2}).copy();

	// Move it to Python IceSheet object.
	self->init(std::move(sheet));

	return 0;
}

static void IceSheet_dealloc(PyIceSheet *self)
{
	self->~PyIceSheet();
	self->ob_type->tp_free((PyObject *)self);
}

//static PyMemberDef IceSheet_members[] = {{NULL}};

// -----------------------------------------------------------

static PyMethodDef IceSheet_methods[] = {
	{NULL}     /* Sentinel - marks the end of this structure */
};

static PyMemberDef IceSheet_members[] = {
	{NULL}
};

PyTypeObject IceSheetType = {
  PyVarObject_HEAD_INIT(NULL, 0)
   "IceSheet",               /* tp_name */
   sizeof(PyIceSheet),     /* tp_basicsize */
   0,                         /* tp_itemsize */
   (destructor)IceSheet_dealloc, /* tp_dealloc */
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
   "IceSheet object",        /* tp_doc */
   0,                         /* tp_traverse */
   0,                         /* tp_clear */
   0,                         /* tp_richcompare */
   0,                         /* tp_weaklistoffset */
   0,                         /* tp_iter */
   0,                         /* tp_iternext */
   IceSheet_methods,         /* tp_methods */
   IceSheet_members,         /* tp_members */
//   0,                         /* tp_members */
   0,                         /* tp_getset */
   0,                         /* tp_base */
   0,                         /* tp_dict */
   0,                         /* tp_descr_get */
   0,                         /* tp_descr_set */
   0,                         /* tp_dictoffset */
   (initproc)IceSheet__init,  /* tp_init */
   0,                         /* tp_alloc */
   (newfunc)&IceSheet_new    /* tp_new */
//   (freefunc)IceSheet_free	/* tp_free */
};
