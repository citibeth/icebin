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

#include <vector>
#include <string>
#include "pyutil.hpp"
#include <cstdio>

namespace giss {

// ===============================================================
PyObject *init_module(
PyModuleDef &module_def,
PyMethodDef **function_sets,	// NULL-terminated Array of NULL-terminated arrays
PyTypeObject **types)			// NULL-terminated Array of type pointers
{
   PyObject* mod;

	// ============ Create flat array of all functions in this module
	int ndefs = 0;
	for (int i=0; function_sets[i] != NULL; ++i) {
		PyMethodDef *functions = function_sets[i];
		for (int j=0; functions[j].ml_name != NULL; ++j) {
			++ndefs;
		}
	}
	PyMethodDef *defs = new PyMethodDef[ndefs+1];	// We will give this memory to Py_InitModule3
	int k=0;
	for (int i=0; function_sets[i] != NULL; ++i) {
		PyMethodDef *functions = function_sets[i];
		for (int j=0; functions[j].ml_name != NULL; ++j) {
			defs[k++] = functions[j];
		}
	}
	const PyMethodDef dummy = {NULL};
	defs[k++] = dummy;


   // ============= Create the module
   module_def.m_methods = defs;

   mod = PyModule_Create(&module_def);
   if (mod == NULL) {
      return NULL;
   }

	// ============= Add the types to the module
	for (int i=0; types[i] != NULL; ++i) {
		PyTypeObject *type = types[i];
		if (PyType_Ready(type) < 0) return NULL;

		// Add the type to the module.
		Py_INCREF(type);
//		printf("Adding type %s (%p) to module %p\n", type->tp_name, type, mod);
		PyModule_AddObject(mod, type->tp_name, (PyObject*)type);
	}

	// ========== Return the module we creaed, for posterity
	return mod;
}


// ===============================================================
/** Makes sure a Numpy Array has the given type and dimension.  Raises
a Python exception if it does not.
@param type_num See http://docs.scipy.org/doc/numpy/reference/c-api.dtype.html eg: NPY_DOUBLE */
PyArrayObject *check_dimensions(
PyObject *ovec,
std::string const &vname,
int type_num,
int ndim,
int const *dims,
int N_ndim)
{
	// Check user and template #-dimensions match
	if (ndim != N_ndim) {
 		char buf[200 + vname.size()];
		sprintf(buf, "check_dimensions(%s): Number of template dimensions (%d) must match number of array dimensions (%d)", vname.c_str(), N_ndim, ndim);
		fprintf(stderr, "%s\n", buf);
		PyErr_SetString(PyExc_ValueError, buf);
		throw std::exception();
	}

	// Check that it's not null
	if (!ovec) {
		std::string serr = "check_dimensions: Array object " + vname + " is null";
		fprintf(stderr, "%s\n", serr.c_str());
		PyErr_SetString(PyExc_ValueError, serr.c_str());
		throw std::exception();
	}

	// Check that it's type PyArrayObject
	if (!PyArray_Check(ovec)) {
		std::string serr = "check_dimensions: Object " + vname + " is not a Numpy array";
		fprintf(stderr, "%s\n", serr.c_str());
		PyErr_SetString(PyExc_ValueError,
			serr.c_str());
		throw std::exception();
	}
	PyArrayObject *vec = (PyArrayObject *)ovec;

	// Check the data type and number of dimensions
	if (vec->descr->type_num != type_num || vec->nd != ndim)  {
		char buf[200 + vname.size()];
		sprintf(buf, "check_dimensions: %s must be of type_num %d and %d dimensions (its is of type_num=%d and %d dimensions).", vname.c_str(), type_num, ndim, vec->descr->type_num, vec->nd);
		fprintf(stderr, "%s\n", buf);
		PyErr_SetString(PyExc_ValueError, buf);
		throw std::exception();
	}

	// Check the dimensions themselves
	for (int i=0; i<ndim; ++i) {
		if (dims[i] < 0) continue;		// Don't check this dimension
		if (dims[i] != vec->dimensions[i]) {
			char buf[200 + vname.size()];
			sprintf(buf,
				"%s: Array dimension #%d is %d, should be %d",
				vname.c_str(), i, vec->dimensions[i], dims[i]);
			fprintf(stderr, "%s\n", buf);
			PyErr_SetString(PyExc_ValueError, buf);
			throw std::exception();
		}
	}
	return vec;
}
// ===============================================================
#if 0
VectorSparseMatrix py_to_VectorSparseMatrix(SparseDescr const &descr, rows_py, cols_py, data_py)
{
	// Check arrays and copy to std::vector
	auto rows(py_to_vector<int>(rows_py, "rows");
	auto cols(py_to_vector<int>(cols_py, "cols", rows.size());
	auto data(py_to_vector<double>(data_py, "data", rows.size());

	return giss::VectorSparseMatrix(descr,
		std::move(rows), std::move(cols), std::move(data));
}
#endif

PyObject *VectorSparseMatrix_to_py(VectorSparseMatrix const &mat)
{
	PyObject *rows_py = NULL;
	PyObject *cols_py = NULL;
	PyObject *data_py = NULL;

	try {
		const int dims[] = {mat.size()};
		const int ndim = 1;

		rows_py = vector_to_py(mat.rows());
		cols_py = vector_to_py(mat.cols());
		data_py = vector_to_py(mat.vals());

		return Py_BuildValue("iiOOO",
			mat.nrow, mat.ncol,
			rows_py, cols_py, data_py);
	} catch(...) {
		// de-allocate since we're not returning it
		if (rows_py) Py_DECREF(rows_py);
		if (cols_py) Py_DECREF(cols_py);
		if (data_py) Py_DECREF(data_py);
		throw std::exception();
	}
}

VectorSparseMatrix py_to_VectorSparseMatrix(PyObject *m_tuple, std::string const &vname)
{
//printf("py_to_VectorSparseMatrix()\n"); fflush(stdout);

	// Get Arguments
	int nrow;
	int ncol;
	PyObject *rows_py;
	PyObject *cols_py;
	PyObject *data_py;
	if (!PyArg_ParseTuple(m_tuple, "iiOOO",
		&nrow, &ncol,
		&rows_py, &cols_py, &data_py))
	{
 		char buf[200 + vname.size()];
		sprintf(buf, "py_to_vectorSparseMatrix(%s): Trouble parsing tuples", vname.c_str());
		fprintf(stderr, "%s\n", buf);
		PyErr_SetString(PyExc_ValueError, buf);

		throw std::exception();
	}
printf("pyutil: nrow=%d, ncol=%d, rows_py=%p, cols_py=%p, data_py=%p\n", nrow, ncol, rows_py, cols_py, data_py);

	// Check arrays and copy to std::vector
	auto rows(py_to_vector<int>(rows_py, "rows"));
	auto cols(py_to_vector<int>(cols_py, "cols", rows.size()));
	auto data(py_to_vector<double>(data_py, "data", rows.size()));

	return giss::VectorSparseMatrix(SparseDescr(nrow, ncol),
		std::move(rows), std::move(cols), std::move(data));
}


giss::BlitzSparseMatrix py_to_BlitzSparseMatrix(PyObject *m_tuple, std::string const &vname)
{
//printf("py_to_BlitzSparseMatrix()\n"); fflush(stdout);

	// Get Arguments
	int nrow;
	int ncol;
	PyObject *rows_py;
	PyObject *cols_py;
	PyObject *data_py;
	if (!PyArg_ParseTuple(m_tuple, "iiOOO",
		&nrow, &ncol,
		&rows_py, &cols_py, &data_py))
	{
 		char buf[200 + vname.size()];
		sprintf(buf, "py_to_BlitzSparseMatrix(%s): Trouble parsing tuples", vname.c_str());
		fprintf(stderr, "%s\n", buf);
		PyErr_SetString(PyExc_ValueError, buf);

		throw std::exception();
	}
//printf("pyutil: nrow=%d, ncol=%d, rows_py=%p, cols_py=%p, data_py=%p\n", nrow, ncol, rows_py, cols_py, data_py);

	// Check arrays and copy to std::vector
	int dims[1] = {-1};
	auto rows(py_to_blitz<int,1>(rows_py, "rows", 1, dims));
	dims[0] = rows.size();
	auto cols(py_to_blitz<int,1>(cols_py, "cols", 1, dims));
	auto data(py_to_blitz<double,1>(data_py, "data", 1, dims));

//printf("py_to_blitz stuff: nnz = %d\n", dims[0]);

	return giss::BlitzSparseMatrix(SparseDescr(nrow, ncol),
		rows, cols, data);
}






}
