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
		char *grid_fname_py;
		char *sproj_py;
		char *sdir_py;
		if (!PyArg_ParseTuple(args, "sss",
			&grid_fname_py, &sproj_py, &sdir_py))
		{ return NULL; }

		// Load the grid from the file
		// (so we don't have to wrap the Grid class in Python for now)
		NcFile nc(grid_fname_py);
		std::unique_ptr<glint2::Grid> grid(glint2::read_grid(nc, "grid"));
		nc.close();

		// Set up the projection
		std::string sproj(sproj_py);
		giss::Proj2 proj;
		if (grid->scoord == "lonlat") {
			proj.init(sproj, giss::Proj2::Direction::LL2XY);
		} else {
			fprintf(stderr, "proj_to_native() only makes sense for grids in Lon/Lat Coordinates!");
			throw std::exception();
		}

		// Do the call
		std::vector<double> ret(proj_native_area_correct(
			*grid, proj, std::string(sdir_py)));

		// Create an output tuple of Numpy arrays
		ret_py = giss::vector_to_py(ret);
		return ret_py;
	} catch(...) {
		if (ret_py) Py_DECREF(ret_py);
		return NULL;
	}
}




PyMethodDef matrix_ops_functions[] = {
	{"coo_matvec", (PyCFunction)coo_matvec_py, METH_VARARGS,
		""},
	{"grid1_to_grid2", (PyCFunction)grid1_to_grid2_py, METH_VARARGS,
		""},
	{"grid2_to_grid1", (PyCFunction)grid1_to_grid2_py, METH_VARARGS,
		""},
	{"mask_out", (PyCFunction)mask_out_py, METH_VARARGS,
		""},
	{"proj_native_area_correct", (PyCFunction)proj_native_area_correct_py, METH_VARARGS,
		""},

	{NULL}     /* Sentinel - marks the end of this structure */
};

