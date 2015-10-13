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

#include <structmember.h>	// Also Python-related

#include <netcdfcpp.h>
#include "PyClass.hpp"

#define RETURN_INVALID_ARGUMENTS(fname) {\
	fprintf(stderr, fname "(): invalid arguments.\n"); \
	PyErr_SetString(PyExc_ValueError, fname "(): invalid arguments."); \
	return NULL; }


// ========================================================================

// ========= class _glint2.NcFile :
/** Read from a file */
static int NcFile__init(PyClass<NcFile> *self, PyObject *args, PyObject *kwds)
{
	try {
		// Get arguments
		char *fname_py = NULL;
		char *filemode_py = NULL;
		if (!PyArg_ParseTuple(args, "ss", &fname_py, &filemode_py)) {
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"Bad arguments to NcFile_save().");
			return 0;
		}

		auto mode = ((strcmp(filemode_py, "w") == 0) ?
			NcFile::Replace : NcFile::ReadOnly);

		// Instantiate C++ NcFile
		std::unique_ptr<NcFile> maker(new NcFile(fname_py, mode));

		// Move it to Python NcFile object.
		self->init(std::move(maker));
		return 0;
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in NcFile__init()");
		return 0;
	}
}


/** Read from a file */
static PyObject *NcFile_close(PyClass<NcFile> *self, PyObject *args)
{
	try {
		// Get arguments
		if (!PyArg_ParseTuple(args, "")) {
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"NcFile_realize() takes no arguments.");
			return 0;
		}
		self->ptr->close();

		return Py_None;
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in NcFile_realize()");
		return 0;
	}
}


// -----------------------------------------------------------

static PyMethodDef NcFile_methods[] = {

	{"close", (PyCFunction)NcFile_close, METH_VARARGS,
		""},
	{NULL}     /* Sentinel - marks the end of this structure */
};

static PyMemberDef NcFile_members[] = {
	{NULL}
};

PyTypeObject NcFileType = {
  PyVarObject_HEAD_INIT(NULL, 0)
   "NcFile",               /* tp_name */
   sizeof(PyClass<NcFile>),     /* tp_basicsize */
   0,                         /* tp_itemsize */
   (destructor)&PyClass<NcFile>::dealloc, /* tp_dealloc */
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
   "NcFile object",        /* tp_doc */
   0,                         /* tp_traverse */
   0,                         /* tp_clear */
   0,                         /* tp_richcompare */
   0,                         /* tp_weaklistoffset */
   0,                         /* tp_iter */
   0,                         /* tp_iternext */
   NcFile_methods,         /* tp_methods */
   NcFile_members,         /* tp_members */
//   0,                         /* tp_members */
   0,                         /* tp_getset */
   0,                         /* tp_base */
   0,                         /* tp_dict */
   0,                         /* tp_descr_get */
   0,                         /* tp_descr_set */
   0,                         /* tp_dictoffset */
   (initproc)NcFile__init,  /* tp_init */
   0,                         /* tp_alloc */
   (newfunc)&PyClass<NcFile>::new_instance    /* tp_new */
//   (freefunc)NcFile_free	/* tp_free */
};
