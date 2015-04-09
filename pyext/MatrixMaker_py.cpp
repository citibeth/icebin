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

#include "MatrixMaker_py.hpp"
#include <structmember.h>	// Also Python-related
#include <numpy/arrayobject.h>
#include <math.h>
#include <glint2/MatrixMaker.hpp>
#include "pyutil.hpp"
#include <giss/memory.hpp>
#include "PyClass.hpp"
#include <glint2/IceSheet_L0.hpp>
#include <giss/enum.hpp>

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

/** Used to construct MatrixMaker */
static int MatrixMaker__init(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{

	try {
		// Get arguments
		int correct_area1 = true;
		static char const *keyword_list[] = {"correct_area1", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "|i",
			const_cast<char **>(keyword_list),
			&correct_area1))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"MatrixMaker__init() called without valid arguments.");
			return 0;
		}

		// Instantiate C++ Ice Maker
		std::unique_ptr<GridDomain> domain(new GridDomain_Identity());
		std::unique_ptr<glint2::MatrixMaker> maker(new MatrixMaker(
			correct_area1, std::move(domain)));

		// Move it to Python MatrixMaker object.
		self->init(std::move(maker));
		return 0;
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker__init()");
		return 0;
	}
}

/** Use this to generate GLINT2 config file */
static PyObject *MatrixMaker_init(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{
	try {
		// Get arguments
		const char *grid1_fname_py = NULL;
		const char *shc_index_type = NULL;
		PyObject *hpdefs_py = NULL;
		PyObject *mask1_py = NULL;
		int correct_area1 = true;

		static char const *keyword_list[] = {"grid1_fname", "hc_index_type", "hpdefs", "mask1", "correct_area1", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "ssO|Oi",
			const_cast<char **>(keyword_list),
			&grid1_fname_py, &shc_index_type, &hpdefs_py,
			&mask1_py, &correct_area1))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"MatrixMaker_init() called without valid arguments.");
			return 0;
		}


		// Instantiate C++ Ice Maker
		std::unique_ptr<GridDomain> domain(new GridDomain_Identity());
		std::unique_ptr<glint2::MatrixMaker> maker(new MatrixMaker(
			correct_area1, std::move(domain)));

		maker->_hc_index_type = giss::parse_enum<HCIndex::Type>(shc_index_type);

		{NcFile nc(grid1_fname_py, NcFile::ReadOnly);
			maker->grid1 = read_grid(nc, "grid");
		}

		int n1 = maker->grid1->ndata();

		if (mask1_py) {
			maker->mask1.reset(new blitz::Array<int,1>(
				giss::py_to_blitz<int,1>(mask1_py, "mask1", {n1})));
		}

		maker->hpdefs = giss::py_to_vector<double>(hpdefs_py, "hpdefs", -1);

		// Move it to Python MatrixMaker object.
		self->init(std::move(maker));
		return Py_None;
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_init()");
		return 0;
	}
}

/** Use this to generate GLINT2 config file */
static PyObject *MatrixMaker_set_mask1(PyMatrixMaker *self, PyObject *args)
{
	try {
		// Get arguments
		PyObject *mask1_py = NULL;

		if (!PyArg_ParseTuple(
			args, "O",
			&mask1_py))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"MatrixMaker_set_mask1() called without valid arguments.");
			return 0;
		}

		self->maker->mask1.reset(new blitz::Array<int,1>(
			giss::py_to_blitz<int,1>(mask1_py, "mask1", {self->maker->n1()})));

		return Py_None;
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_set_mask1()");
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
		const char *name_py = NULL;

		static char const *keyword_list[] = {"grid2_fname", "exgrid_fname", "elev2", "mask2", "name", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "ssO|Os",
			const_cast<char **>(keyword_list),
			&grid2_fname_py, &exgrid_fname_py,
			&elev2_py, &mask2_py, &name_py))
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
		if (name_py) sheet->name = name_py;
		sheet->grid2 = std::move(grid2);
		int n2 = sheet->grid2->ndata();

		// Downcast to ExchangeGrid
		{NcFile nc(exgrid_fname_py, NcFile::ReadOnly);
		sheet->exgrid = giss::dynamic_cast_unique_ptr<ExchangeGrid>(
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

static PyObject *MatrixMaker_hp_to_iceinterp(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{
	PyObject *ret_py = NULL;
	try {
		// Get arguments
		const char *ice_sheet_name_py;
		const char *dest_py = "ICE";
		int fill_masked = 0;
		static char const *keyword_list[] = {"sheetname", "dest", "fill_masked", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "s|si",
			const_cast<char **>(keyword_list),
			&ice_sheet_name_py, &dest_py, &fill_masked))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"hp_to_iceinterp() called without a valid string as argument.");
			return 0;
		}
		glint2::MatrixMaker *maker = self->maker.get();
		std::string const ice_sheet_name(ice_sheet_name_py);
		auto dest(giss::parse_enum<IceInterp>(dest_py));

		// Look up the ice sheet
		IceSheet *sheet = (*maker)[ice_sheet_name];
		if (!sheet) {
			PyErr_SetString(PyExc_ValueError,
				("Could not find ice sheet named " + ice_sheet_name).c_str());
			return 0;
		}

		// Get the hp_to_iceinterp matrix from it
		auto ret_c(sheet->hp_to_iceinterp(dest, fill_masked));

		// Create an output tuple of Numpy arrays
		ret_py = giss::VectorSparseMatrix_to_py(*ret_c);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_hp_to_iceinterp()");
		return 0;
	}
}

static PyObject *MatrixMaker_iceinterp_to_atm(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{
	PyObject *ret_py = NULL;
	try {

		// Get arguments
		char *ice_sheet_name_py;
		char *src_py = "ICE";
		static char const *keyword_list[] = {"sheetname", "src", NULL};


		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "s|s",
			const_cast<char **>(keyword_list),
			&ice_sheet_name_py, &src_py))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"ice_to_atm() called without a valid string as argument.");
			return 0;
		}
		glint2::MatrixMaker *maker = self->maker.get();
		std::string const ice_sheet_name(ice_sheet_name_py);
		auto src(giss::parse_enum<IceInterp>(src_py));

		// Look up the ice sheet
		IceSheet *sheet = (*maker)[ice_sheet_name];
		if (!sheet) {
			PyErr_SetString(PyExc_ValueError,
				("Could not find ice sheet named " + ice_sheet_name).c_str());
			return 0;
		}

		// Get the ice_to_atm matrix from it
		giss::MapSparseVector<int,double> area1_m, area1_m_inv;
		auto ret_c(sheet->iceinterp_to_projatm(area1_m, src));
		if (maker->correct_area1)
			sheet->atm_proj_correct(area1_m, ProjCorrect::PROJ_TO_NATIVE);
		divide_by(*ret_c, area1_m, area1_m_inv);

		// Create an output tuple of Numpy arrays
		ret_py = giss::VectorSparseMatrix_to_py(*ret_c);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_ice_to_atm()");
		return 0;
	}
}



static PyObject *MatrixMaker_hp_to_atm(PyMatrixMaker *self, PyObject *args)
{
	PyObject *ret_py = NULL;
	try {

		// Get arguments
		if (!PyArg_ParseTuple(args, "")) {
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"hp_to_atm() does not take arguments.");
			return 0;
		}

		// Get the hp_to_atm matrix from it
		auto ret_c(self->maker->hp_to_atm());

		// Create an output tuple of Numpy arrays
		ret_py = giss::VectorSparseMatrix_to_py(*ret_c);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_hp_to_atm()");
		return 0;
	}
}

static PyObject *MatrixMaker_atm_to_hp(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{
	PyObject *ret_py = NULL;
	try {

		PyObject *f1_py;
		int force_lambda = 0;	// false
		static char const *keyword_list[] = {"f1", "force_lambda", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "O|i",
			const_cast<char **>(keyword_list),
			&f1_py, &force_lambda))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"atm_to_hp(): bad arguments");
			return 0;
		}

		int dims[1] = {self->maker->n1()};
		auto f1(giss::py_to_blitz<double,1>(f1_py, "f1", 1, dims));

		// Call!
		giss::VectorSparseVector<int, double> f3(self->maker->atm_to_hp(f1, force_lambda));

		// Copy output for return
		blitz::Array<double,1> ret(self->maker->n3());
		to_blitz(f3, ret);
#if 0
		ret = 0;
		for (auto ii = f3.begin(); ii != f3.end(); ++ii) {
			int i3 = ii->first;
			double val = ii->second;
			ret(i3) = val;
		}
#endif

		ret_py = giss::blitz_to_py(ret);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_hp_to_atm()");
		return 0;
	}
}

static PyObject *MatrixMaker_iceinterp_to_hp(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{
printf("BEGIN MatrixMaker_ice_to_hp()\n");
	PyObject *ret_py = NULL;
	try {

		// Get arguments
		PyObject *f2s_py;		// Should be [(1 : [...]), (2 : [...])]
		PyObject *initial_py = NULL;		// array[n3]
		const char *src_py = "ICE";
		const char *qp_algorithm_py = "SINGLE_QP";
		static char const *keyword_list[] = {"f2s", "initial3", "src", "qp_algorithm", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "O|Oss",
			const_cast<char **>(keyword_list),
			&f2s_py,
			&initial_py, &src_py, &qp_algorithm_py)) {
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"Bad arguments for ice_to_hp().");
			return 0;
		}

		auto src(giss::parse_enum<IceInterp>(src_py));
		auto qp_algorithm(giss::parse_enum<QPAlgorithm>(qp_algorithm_py));

		if (!PyList_Check(f2s_py)) {
			PyErr_SetString(PyExc_ValueError,
				"Argument must be a list.");
			return 0;
		}


		// Convert the Python dict to a C++ dict
		std::map<int, blitz::Array<double,1>> f2s;
		Py_ssize_t len = PyList_Size(f2s_py);
		for (int i=0; i < len; ++i) {
			PyObject *ii = PyList_GetItem(f2s_py, i);
			if (!PyTuple_Check(ii)) {
				PyErr_SetString(PyExc_ValueError,
					"List must contain tuples");
				return 0;
			}
			char *sheetname_py;
			PyObject *f2_py;
			PyArg_ParseTuple(ii, "sO", &sheetname_py, &f2_py);
printf("MatrixMaker_ice_to_hp(): Adding %s\n", sheetname_py);
			IceSheet *sheet = (*self->maker)[std::string(sheetname_py)];

			int dims[1] = {src == IceInterp::ICE ? sheet->n2() : sheet->n4()};
			auto f2(giss::py_to_blitz<double,1>(f2_py, "f2", 1, dims));

			f2s.insert(std::make_pair(sheet->index, f2));

		}

		// Get the blitz array from python
		int dims[1] = {self->maker->n3()};
		blitz::Array<double,1> initial;
		if (initial_py) {
			initial.reference(giss::py_to_blitz<double,1>(initial_py, "initial", 1, dims));
		}

		// Call!
		giss::VectorSparseVector<int, double> f3(
			self->maker->iceinterp_to_hp(f2s, initial, src, qp_algorithm));

		// Copy output for return
		blitz::Array<double,1> ret(self->maker->n3());
		giss::to_blitz(f3, ret);
#if 0
		ret = 0;
		for (auto ii = f3.begin(); ii != f3.end(); ++ii) {
			int i3 = ii->first;
			double val = ii->second;
			ret(i3) = val;
		}
#endif

		ret_py = giss::blitz_to_py(ret);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_ice_to_hp()");
		return 0;
	}
}

static PyObject *MatrixMaker_ice_to_interp(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{
printf("BEGIN MatrixMaker_ice_to_interp()\n");
	PyObject *ret_py = NULL;
	try {

		// Get arguments
		const char *ice_sheet_name_py;
		PyObject *f2_py;		// Should be [(1 : [...]), (2 : [...])]
		static char const *keyword_list[] = {"sheetname", "f2", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "sO",
			const_cast<char **>(keyword_list),
			&ice_sheet_name_py, &f2_py))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"Bad arguments for ice_to_interp().");
			return 0;
		}

		glint2::MatrixMaker *maker = self->maker.get();
		std::string const ice_sheet_name(ice_sheet_name_py);
		IceSheet *sheet = (*maker)[ice_sheet_name];

		int dims[1] = {sheet->n2()};
		auto f2(giss::py_to_blitz<double,1>(f2_py, "f2", 1, dims));

		auto ret(sheet->ice_to_interp(f2));

		ret_py = giss::blitz_to_py(ret);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_ice_to_hp()");
		return 0;
	}
}


static PyObject *MatrixMaker_area1(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{
printf("BEGIN MatrixMaker_area1()\n");
	PyObject *ret_py = NULL;
	try {

		// Get arguments
		const char *ice_sheet_name_py;
		static char const *keyword_list[] = {"sheetname", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "s",
			const_cast<char **>(keyword_list),
			&ice_sheet_name_py))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"Bad arguments for area1().");
			return 0;
		}

		glint2::MatrixMaker *maker = self->maker.get();
		std::string const ice_sheet_name(ice_sheet_name_py);
		IceSheet *sheet = (*maker)[ice_sheet_name];

		giss::MapSparseVector<int,double> area1_m;
		sheet->accum_areas(area1_m);

		blitz::Array<double,1> ret(maker->n1());
		ret = 0;
		for (auto ii=area1_m.begin(); ii != area1_m.end(); ++ii) {
			ret(ii->first) = ii->second;
		}


		ret_py = giss::blitz_to_py(ret);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_ice_to_hp()");
		return 0;
	}
}

static PyObject *MatrixMaker_set_interp_grid(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{
	PyObject *ret_py = NULL;
	try {
		// Get arguments
		const char *interp_grid_py = "EXCH";
		static char const *keyword_list[] = {"interp_grid", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "s",
			const_cast<char **>(keyword_list),
			&interp_grid_py))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"set_interp_grid() called without a valid string as argument.");
			return 0;
		}
		glint2::MatrixMaker *maker = self->maker.get();
		auto interp_grid(giss::parse_enum<IceExch>(interp_grid_py));

		for (auto sheet = maker->sheets.begin(); sheet != maker->sheets.end(); ++sheet) {
			IceSheet_L0 *sheet0 = dynamic_cast<IceSheet_L0 *>(&*sheet);
			if (sheet0) sheet0->interp_grid = interp_grid;
		}

		return Py_None;
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_hp_to_iceinterp()");
		return 0;
	}
}

// ---------------------------------------------------------
static PyObject *MatrixMaker_set_interp_style(PyMatrixMaker *self, PyObject *args, PyObject *kwds)
{
	PyObject *ret_py = NULL;
	try {
		// Get arguments
		const char *interp_style_py = "Z_INTERP";
		static char const *keyword_list[] = {"interp_style", NULL};

		if (!PyArg_ParseTupleAndKeywords(
			args, kwds, "s",
			const_cast<char **>(keyword_list),
			&interp_style_py))
		{
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"set_interp_style() called without a valid string as argument.");
			return 0;
		}
		glint2::MatrixMaker *maker = self->maker.get();
		auto interp_style(giss::parse_enum<InterpStyle>(interp_style_py));

		for (auto sheet = maker->sheets.begin(); sheet != maker->sheets.end(); ++sheet) {
			sheet->interp_style = interp_style;
		}

		return Py_None;
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_hp_to_iceinterp()");
		return 0;
	}
}

// ---------------------------------------------------------


/** Read from a file */
static PyObject *MatrixMaker_write(PyMatrixMaker *self, PyObject *args)
{
	try {
		// Get arguments
		PyClass<NcFile> *nc_py = NULL;
		char *vname_py = NULL;
		if (!PyArg_ParseTuple(args, "Os", &nc_py, &vname_py)) {
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"Bad arguments to MatrixMaker_write().");
			return 0;
		}

//		NcFile nc(fname_py, NcFile::Replace);
		self->maker->netcdf_define(*(nc_py->ptr), std::string(vname_py))();
//		nc.close();

		return Py_None;
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_write()");
		return 0;
	}
}

/** Read from a file */
static PyObject *MatrixMaker_load(PyMatrixMaker *self, PyObject *args)
{
	try {
		// Get arguments
		char *fname_py = NULL;
		char *vname_py = NULL;
		if (!PyArg_ParseTuple(args, "ss", &fname_py, &vname_py)) {
			// Throw an exception...
			PyErr_SetString(PyExc_ValueError,
				"Bad arguments to MatrixMaker_load().");
			return 0;
		}

		NcFile nc(fname_py, NcFile::ReadOnly);
		self->maker->read_from_netcdf(nc, std::string(vname_py));
		nc.close();

		return Py_None;
	} catch(...) {
		PyErr_SetString(PyExc_ValueError, "Error in MatrixMaker_load()");
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

	{"init", (PyCFunction)MatrixMaker_init, METH_KEYWORDS,
		"Function used to contruct a MatrixMaker from scratch\n"
		"(rather than loading from a file.)\n"
		"\n"
		"grid1_fname : str\n"
		"    Name of the grid file containing the GCM grid\n"
		"\n"
		"hc_index_type : str (see HCIndex::Type)\n"
		"    Indexing convention for elevation class arrays where:\n"
		"        nhp = index for elevation class (zero-based)\n"
		"        n1 = 1D collapsed index of GCM grid cell (zero-based)\n"
		"    'MODELE' : A[nhp, n1]\n"
		"\n"
		"hpdefs : double[nhp] (m)\n"
		"    The elevation points to use.  The same elevation points will be\n"
		"    used for all GCM grid cells.\n"
		"\n"
		"mask1 : int[n1] (OPTIONAL)\n"
		"    Use to mark GCM cells as not used.  Masking convention is the same\n"
		"    as in Matplotlib: zero means the cell is OK, non-zero means it is\n"
		"    not used.\n"
		"\n"
		"correct_area1 : int/bool (DEFAULT 1)\n"
		"    if non-zero: account for projection error in transformations.\n"
	}
,

	{"add_ice_sheet", (PyCFunction)MatrixMaker_add_ice_sheet, METH_KEYWORDS,
		""},
	{"hp_to_iceinterp", (PyCFunction)MatrixMaker_hp_to_iceinterp, METH_KEYWORDS,
		""},
	{"hp_to_atm", (PyCFunction)MatrixMaker_hp_to_atm, METH_KEYWORDS,
		""},
	{"atm_to_hp", (PyCFunction)MatrixMaker_atm_to_hp, METH_KEYWORDS,
		""},
	{"iceinterp_to_atm", (PyCFunction)MatrixMaker_iceinterp_to_atm, METH_KEYWORDS,
		""},
	{"iceinterp_to_hp", (PyCFunction)MatrixMaker_iceinterp_to_hp, METH_KEYWORDS,
		""},
	{"ice_to_interp", (PyCFunction)MatrixMaker_ice_to_interp, METH_KEYWORDS,
		""},
	{"area1",  (PyCFunction)MatrixMaker_area1, METH_KEYWORDS,
		""},
	{"set_interp_grid",  (PyCFunction)MatrixMaker_set_interp_grid, METH_KEYWORDS,
		""},
	{"set_interp_style",  (PyCFunction)MatrixMaker_set_interp_style, METH_KEYWORDS,
		""},
	{"set_mask1",  (PyCFunction)MatrixMaker_set_mask1, METH_VARARGS,
		""},
	{"realize", (PyCFunction)MatrixMaker_realize, METH_VARARGS,
		"Call this after constructing a MatrixMaker with init() and\n"
		"add_ice_sheet() methods.  Takes no arguments are returns\n"
		"nothing."},
	{"write", (PyCFunction)MatrixMaker_write, METH_VARARGS,
		"Write the MatrixMaker structure to a netCDF file."
		"\n"
		"nc : netCDF4.Dataset\n"
		"    The open netCDF file\n"
		"vname : str\n"
		"    The name of the variable to write it to.\n"},
	{"load", (PyCFunction)MatrixMaker_load, METH_VARARGS,
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
