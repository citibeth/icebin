#define NO_IMPORT_ARRAY
#include "_glint2_module.hpp"

#include "pyutil.hpp"
#include <glint2/matrix_ops.hpp>
#include <giss/SparseMatrix.hpp>
#include "Grid_py.hpp"

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

PyObject *height_classify_py(PyObject *self, PyObject *args)
{
	PyObject *ret_py = NULL;
	try {
		// Get Arguments
		PyObject *overlap_py;
		PyObject *elev2_py;
		PyObject *hcmax_py;
		if (!PyArg_ParseTuple(args, "OOO",
			&overlap_py, &elev2_py, &hcmax_py))
		{ return NULL; }

		// Check arrays and copy to giss::VectorSparseMatrix
		auto overlap(giss::py_to_BlitzSparseMatrix(overlap_py, "overlap"));
		int dims[1] = {overlap.ncol};	// |grid2|
		auto elev2(giss::py_to_blitz<double,1>(elev2_py, "elev2", 1, dims));
		dims[0] = -1;
		auto hcmax(giss::py_to_blitz<double,1>(hcmax_py, "hcmax", 1, dims));

		// Do the call
		std::unique_ptr<giss::VectorSparseMatrix> ret_c(
			glint2::height_classify(overlap, elev2, hcmax));

		// Create an output tuple of Numpy arrays
		ret_py = giss::VectorSparseMatrix_to_py(*ret_c);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		return NULL;
	}
}

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

		// Keep track of which items we've written to.
		std::vector<bool> written(mat.nrow, false);

		// Do the multiplication, and we're done!
		int nnz = mat.size();
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
		auto overlap(giss::py_to_BlitzSparseMatrix(overlap_py, "overlap"));

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

PyObject *grid2_to_grid1_py(PyObject *self, PyObject *args)
{
	PyObject *ret_py = NULL;
	try {
		// Get Arguments
		PyObject *overlap_py;
		if (!PyArg_ParseTuple(args, "O",
			&overlap_py))
		{ return NULL; }

		// Check arrays and copy to giss::VectorSparseMatrix
		auto overlap(giss::py_to_BlitzSparseMatrix(overlap_py, "overlap"));

		// Do the call
		std::unique_ptr<giss::VectorSparseMatrix> ret_c(
			glint2::grid2_to_grid1(overlap));

		// Create an output tuple of Numpy arrays
		ret_py = giss::VectorSparseMatrix_to_py(*ret_c);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		return NULL;
	}
}

PyObject *mask_out_py(PyObject *self, PyObject *args)
{
	PyObject *ret_py = NULL;
	try {
		// Get Arguments
		PyObject *overlap_py;
		PyObject *mask1_py;
		PyObject *mask2_py;
		if (!PyArg_ParseTuple(args, "OOO",
			&overlap_py, &mask1_py, &mask2_py))
		{ return NULL; }

		// Check arrays and copy to giss::VectorSparseMatrix
		auto overlap(giss::py_to_BlitzSparseMatrix(overlap_py, "overlap"));

		blitz::Array<int,1> mask1, mask2;

		printf("mask1_py = %p\n", mask1_py);
		printf("mask2_py = %p\n", mask2_py);

		blitz::Array<int,1> *mask1p = NULL;
		if (mask1_py != Py_None) {
			mask1.reference(giss::py_to_blitz<int,1>(mask1_py, "mask1", {overlap.nrow}));
			mask1p = &mask1;
		}

		blitz::Array<int,1> *mask2p = NULL;
		if (mask2_py != Py_None) {
			mask2.reference(giss::py_to_blitz<int,1>(mask2_py, "mask2", {overlap.ncol}));
			mask2p = &mask2;
		}

		// Do the call
printf("Calling mask_out()\n");
		std::unique_ptr<giss::VectorSparseMatrix> ret_c(
			glint2::mask_out(overlap, mask1p, mask2p));
printf("Done Calling mask_out()\n");

		// Create an output tuple of Numpy arrays
		ret_py = giss::VectorSparseMatrix_to_py(*ret_c);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		return NULL;
	}
}


PyObject *proj_native_area_correct_py(PyObject *self, PyObject *args)
{
	PyObject *ret_py = NULL;
	try {
		// Get Arguments
		PyGrid *grid_py;
		char *sproj_py;
		char *sdir_py;
		if (!PyArg_ParseTuple(args, "Oss",
			&grid_py, &sproj_py, &sdir_py))
		{ return NULL; }

		glint2::Grid &grid(*grid_py->grid);

		// Do the call
		std::vector<double> ret(proj_native_area_correct(
			grid, std::string(sproj_py), std::string(sdir_py)));

		// Create an output tuple of Numpy arrays
		ret_py = giss::vector_to_py(ret);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		return NULL;
	}
}

PyObject *multiply_bydiag_py(PyObject *self, PyObject *args)
{
	PyObject *ret_py = NULL;
	try {
		// Get Arguments
		PyObject *arg1_py;
		PyObject *arg2_py;
		if (!PyArg_ParseTuple(args, "OO",
			&arg1_py, &arg2_py))
		{ return NULL; }
printf("arg1=%p, arg2=%p\n", arg1_py, arg2_py);

		// Figure whether we're doing (A * diag) or (diag * A)
		PyObject *mat_py;
		PyObject *diag_py;
		bool right_mul = PyTuple_Check(arg1_py);
printf("right_mul = %d\n", right_mul);
		if (right_mul) {	// A * diag
			mat_py = arg1_py;
			diag_py = arg2_py;
		} else {
			mat_py = arg2_py;
			diag_py = arg1_py;
		}
printf("mat=%p, diag=%p\n", arg1_py, arg2_py);

		// Check arrays and copy to giss::VectorSparseMatrix
		auto mat(giss::py_to_BlitzSparseMatrix(mat_py, "mat"));

		int dims[1] = {right_mul ? mat.ncol : mat.nrow};
		auto diag(giss::py_to_blitz<double,1>(diag_py, "diag", 1, dims));

		// Do the call
		if (right_mul) multiply_bydiag(mat, diag);
		else multiply_bydiag(diag, mat);

		return Py_None;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		return NULL;
	}
}

PyObject *hp_interp_py(PyObject *self, PyObject *args)
{
	PyObject *ret_py = NULL;
	try {
		// Get Arguments
		PyObject *overlap_py;
		PyObject *elev2_py;
		PyObject *hpdefs_py;
		if (!PyArg_ParseTuple(args, "OOO",
			&overlap_py, &elev2_py, &hpdefs_py))
		{ return NULL; }

		// Check arrays and copy to giss::VectorSparseMatrix
		auto overlap(giss::py_to_BlitzSparseMatrix(overlap_py, "overlap"));
		int dims[1] = {overlap.ncol};	// |grid2|
		auto elev2(giss::py_to_blitz<double,1>(elev2_py, "elev2", 1, dims));
		auto hpdefs(giss::py_to_vector<double>(hpdefs_py, "hpdefs", -1));

		// Do the call
		std::unique_ptr<giss::VectorSparseMatrix> ret_c(
			glint2::hp_interp(overlap, elev2, hpdefs));

		// Create an output tuple of Numpy arrays
		ret_py = giss::VectorSparseMatrix_to_py(*ret_c);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		return NULL;
	}
}




PyMethodDef matrix_ops_functions[] = {
	{"height_classify", (PyCFunction)height_classify_py, METH_VARARGS,
		""},
	{"coo_matvec", (PyCFunction)coo_matvec_py, METH_VARARGS,
		"Compute M*x, taking care with unspecified elements in M"},
	{"grid1_to_grid2", (PyCFunction)grid1_to_grid2_py, METH_VARARGS,
		""},
	{"grid2_to_grid1", (PyCFunction)grid2_to_grid1_py, METH_VARARGS,
		""},
	{"mask_out", (PyCFunction)mask_out_py, METH_VARARGS,
		""},
	{"proj_native_area_correct", (PyCFunction)proj_native_area_correct_py, METH_VARARGS,
		""},
	{"multiply_bydiag", (PyCFunction)multiply_bydiag_py, METH_VARARGS,
		""},
	{"hp_interp", (PyCFunction)hp_interp_py, METH_VARARGS,
		""},

	{NULL}     /* Sentinel - marks the end of this structure */
};

