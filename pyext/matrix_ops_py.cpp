#define NO_IMPORT_ARRAY
#include "_glint2_module.hpp"

#include "pyutil.hpp"
#include <glint2/matrix_ops.hpp>

/*
 * http://dsnra.jpl.nasa.gov/software/Python/numpydoc/numpy-14.html
 * 
 * 
 *  // Make a new double vector of same dimension
 *  101     vecout=(PyArrayObject *) PyArray_FromDims(1,dims,NPY_DOUBLE);
 * 
 * 
 * return PyArray_Return(pyarray)
 * 
 * 
 * return Py_BuildValue("OOO", pyarray1, pyarray2, pyarray3)
 */

PyObject *grid1_to_grid2_py(PyObject *self, PyObject *args)
{
	PyObject *ret_py = NULL;
	try {
		// Get Arguments
		PyObject *overlap_py;
		if (!PyArg_ParseTuple(args, "O",
			&overlap_py))
		{ return NULL; }

		// Check arrays and copy to giss::VectorSparseMatrix
		giss::VectorSparseMatrix *ppp;

		auto overlap(giss::py_to_VectorSparseMatrix(overlap_py, "overlap"));

		// Do the call
		std::unique_ptr<giss::VectorSparseMatrix> ret_c(
			glint2::grid1_to_grid2(overlap));

		// Create an output tuple of Numpy arrays
		ret_py = giss::VectorSparseMatrix_to_py(*ret_c);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		return NULL;
	}
}




PyMethodDef matrix_ops_functions[] = {
	{"grid1_to_grid2", (PyCFunction)grid1_to_grid2_py, METH_VARARGS,
		""},

	{NULL}     /* Sentinel - marks the end of this structure */
};

