#define NO_IMPORT_ARRAY
#include "_glint2_module.hpp"

#include "MatrixMaker_py.hpp"
#include <structmember.h>	// Also Python-related
#include <numpy/arrayobject.h>
#include <math.h>
#include <glint2/MatrixMaker.hpp>
#include "pyutil.hpp"
#include <giss/memory.hpp>

using namespace glint2;

#define RETURN_INVALID_ARGUMENTS(fname) {\
	fprintf(stderr, fname "(): invalid arguments.\n"); \
	PyErr_SetString(PyExc_ValueError, fname "(): invalid arguments."); \
	return NULL; }


// ========================================================================

void PyMatrixMaker::init(std::unique_ptr<glint2::MatrixMaker> &&_maker)
{
	maker = std::move(_maker);
}

// ========= class _glint2.MatrixMaker :
static PyObject *MatrixMaker_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	PyMatrixMaker *self;

	self = new (type->tp_alloc(type, 0)) PyMatrixMaker;

//	Py_INCREF(self);
    return (PyObject *)self;
}

/** Read from a file */
static int MatrixMaker__init(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{
	try {
		// Get arguments
		const char *grid1_fname_py = NULL;
		PyObject *hpdefs_py = NULL;
		PyObject *hcmax_py = NULL;
		PyObject *mask1_py = NULL;

		static char const *keyword_list[] = {"grid1_fname", "hpdefs", "hcmax", "mask1", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "sOO|O",
			const_cast<char **>(keyword_list),
			&grid1_fname_py, &hpdefs_py, &hcmax_py, &mask1_py))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"MatrixMaker__init() called without valid arguments.");
			return 0;
		}


		// Instantiate C++ Ice Maker
		std::unique_ptr<glint2::MatrixMaker> maker(new MatrixMaker);

		{NcFile nc(grid1_fname_py, NcFile::ReadOnly);
			maker->grid1 = read_grid(nc, "grid");
		}

		int n1 = maker->grid1->ndata();

		if (mask1_py) maker->mask1.reset(new blitz::Array<int,1>(
			giss::py_to_blitz<int,1>(mask1_py, "mask1", {n1})));

		maker->hpdefs = giss::py_to_vector<double>(hpdefs_py, "hpdefs", -1);
		int nhc = maker->hpdefs.size();

		maker->hcmax.reference(giss::py_to_blitz<double,1>(hcmax_py, "hcmax", {nhc-1}));

		// Move it to Python MatrixMaker object.
		self->init(std::move(maker));
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker__init()");
		return 0;
	}
}

static PyObject *MatrixMaker_add_ice_sheet(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{
	try {
		// Get arguments
		const char *grid2_fname_py;
		const char *exgrid_fname_py;
		PyObject *elev2_py = NULL;
		PyObject *mask2_py = NULL;

		static char const *keyword_list[] = {"grid2_fname", "exgrid_fname", "elev2", "mask2", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "ssO|O",
			const_cast<char **>(keyword_list),
			&grid2_fname_py, &exgrid_fname_py,
			&elev2_py, &mask2_py))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"MatrixMaker_add_ice_sheet() called without valid arguments.");
			return 0;
		}

		// ========================== Create an IceSheet
		NcFile nc(grid2_fname_py, NcFile::ReadOnly);
		std::unique_ptr<Grid> grid2(read_grid(nc, "grid"));
		nc.close();

		// Instantiate C++ Ice Sheet
		std::unique_ptr<glint2::IceSheet> sheet(
			new_ice_sheet(grid2->parameterization));

		// Fill it in...
		sheet->grid2 = std::move(grid2);
		int n2 = sheet->grid2->ndata();

		// Cast up to ExchangeGrid
		{NcFile nc(exgrid_fname_py, NcFile::ReadOnly);
		sheet->exgrid = giss::unique_cast<ExchangeGrid, Grid>(
			read_grid(nc, "grid"));
		}

		if (mask2_py) sheet->mask2.reset(new blitz::Array<int,1>(
			giss::py_to_blitz<int,1>(mask2_py, "mask2", {n2})));

		sheet->elev2.reference(giss::py_to_blitz<double,1>(elev2_py, "elev2", {n2}));
		// ====================================================

		int ice_sheet_num = self->maker->add_ice_sheet(std::move(sheet));
		return PyInt_FromLong(ice_sheet_num);
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_add_ice_sheet()");
		return 0;
	}
}

/** Read from a file */
static PyObject *MatrixMaker_realize(PyMatrixMaker *self, PyObject *args)
{
	try {
		// Get arguments
		if (!PyArg_ParseTuple(args, "")) {
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"MatrixMaker_realize() takes no arguments.");
			return 0;
		}
		self->maker->realize();

		return Py_None;
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_realize()");
		return 0;
	}
}

static void MatrixMaker_dealloc(PyMatrixMaker *self)
{
	self->~PyMatrixMaker();
	self->ob_type->tp_free((PyObject *)self);
}

//static PyMemberDef MatrixMaker_members[] = {{NULL}};

// -----------------------------------------------------------

static PyMethodDef MatrixMaker_methods[] = {

	{"add_ice_sheet", (PyCFunction)MatrixMaker_add_ice_sheet, METH_KEYWORDS,
		""},
	{"realize", (PyCFunction)MatrixMaker_realize, METH_VARARGS,
		""},
	{NULL}     /* Sentinel - marks the end of this structure */
};

static PyMemberDef MatrixMaker_members[] = {
	{NULL}
};

PyTypeObject MatrixMakerType = {
   PyObject_HEAD_INIT(NULL)
   0,                         /* ob_size */
   "MatrixMaker",               /* tp_name */
   sizeof(PyMatrixMaker),     /* tp_basicsize */
   0,                         /* tp_itemsize */
   (destructor)MatrixMaker_dealloc, /* tp_dealloc */
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
   "MatrixMaker object",        /* tp_doc */
   0,                         /* tp_traverse */
   0,                         /* tp_clear */
   0,                         /* tp_richcompare */
   0,                         /* tp_weaklistoffset */
   0,                         /* tp_iter */
   0,                         /* tp_iternext */
   MatrixMaker_methods,         /* tp_methods */
   MatrixMaker_members,         /* tp_members */
//   0,                         /* tp_members */
   0,                         /* tp_getset */
   0,                         /* tp_base */
   0,                         /* tp_dict */
   0,                         /* tp_descr_get */
   0,                         /* tp_descr_set */
   0,                         /* tp_dictoffset */
   (initproc)MatrixMaker__init,  /* tp_init */
   0,                         /* tp_alloc */
   (newfunc)&MatrixMaker_new    /* tp_new */
//   (freefunc)MatrixMaker_free	/* tp_free */
};
