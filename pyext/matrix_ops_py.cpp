#define NO_IMPORT_ARRAY
#include "_glint2_module.hpp"

#include "pyutil.hpp"
#include <glint2/matrix_ops.hpp>
#include <giss/SparseMatrix.hpp>

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

//static double const nan = std::numeric_limits<double>::quiet_NaN();

PyObject *coo_matvec_py(PyObject *self, PyObject *args)
{
	try {
		PyObject *mat_py = NULL;
		PyObject *xx_py = NULL;
		PyObject *yy_py = NULL;
		if (!PyArg_ParseTuple(args, "OOO",
			&mat_py, &xx_py, &yy_py))
		{ return NULL; }

		// Cast and typecheck arguments
		giss::BlitzSparseMatrix mat(giss::py_to_BlitzSparseMatrix(mat_py, "mat"));
		auto xx(giss::py_to_blitz<double,1>(xx_py, "xx", {mat.ncol}));
		auto yy(giss::py_to_blitz<double,1>(yy_py, "yy", {mat.nrow}));
printf("after py_to_blitz: &yy(0) = %p\n", &yy(0));

		// Keep track of which items we've written to.
		std::vector<bool> written(mat.nrow, false);

		// Do the multiplication, and we're done!
		int nnz = mat.size();
printf("nnz = %d\n", nnz);
		for (int n=0; n<nnz; ++n) {
			int row = mat.rows()(n);
			int col = mat.cols()(n);
			double val = mat.vals()(n);

			// Just do Snowdrift-style "REPLACE".  "MERGE" was never used.
			double old_yy;
			if (written[row]) {
				old_yy = yy(row);
			} else {
				old_yy = 0;
				written[row] = true;
			}
			yy(row) = old_yy + val * xx(col);
printf("%d %d %f (%f)\n", row, col, val, yy(row));
		}
	} catch(...) {
		return NULL;	// Error
	}
	return Py_BuildValue("i", 0);

}


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
	{"coo_matvec", (PyCFunction)coo_matvec_py, METH_VARARGS,
		""},

	{NULL}     /* Sentinel - marks the end of this structure */
};

