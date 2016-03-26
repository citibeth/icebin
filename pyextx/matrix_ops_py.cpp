/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
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


PyObject *coo_matvec_py(PyObject *self, PyObject *args, PyObject *kwds)
{
    try {
        PyObject *mat_py = NULL;
        PyObject *xx_py = NULL;
        PyObject *yy_py = NULL;
        int ignore_nan = 0;

        static char const *keyword_list[] =
            {"mat", "xx", "yy", "ignore_nan", NULL};
        if (!PyArg_ParseTupleAndKeywords(
            args, kwds, "OOO|i",
            const_cast<char **>(keyword_list),
            &mat_py, &xx_py, &yy_py, &ignore_nan))
        {
            PyErr_SetString(PyExc_ValueError,
                "coo_matvec_py() called with invalid arguments.");
            return NULL;
        }

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

            // Ignore NaN in input vector
            if ((ignore_nan != 0) && std::isnan(xx(col))) {
                continue;
            }

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
        return NULL;    // Error
    }
    return Py_BuildValue("i", 0);

}




PyMethodDef matrix_ops_functions[] = {
    {"coo_matvec", (PyCFunction)coo_matvec_py, METH_VARARGS|METH_KEYWORDS,
        "Compute M*x, taking care with unspecified elements in M"},

    {NULL}     /* Sentinel - marks the end of this structure */
};

