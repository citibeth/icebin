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

#pragma once

#include <cassert>
#include <vector>
#include <exception>
#include <blitz/array.h>
#include <giss/SparseMatrix.hpp>
#include <giss/exit.hpp>

namespace giss {

PyObject *init_module(
PyModuleDef &module_def,
PyMethodDef **function_sets,    // NULL-terminated Array of NULL-terminated arrays
PyTypeObject **types);          // NULL-terminated Array of type pointers

// ====================================================================
// --------------------------------------------------------------------
// Convert template types to Numpy type_nums

template<class T>
inline int get_type_num()
{
    PyErr_SetString(PyExc_ValueError, "get_type_num(): Unknown type_num");
    giss::exit(1);
}

template<> inline int get_type_num<short>()
    { return NPY_SHORT; }

template<> inline int get_type_num<int>()
    { return NPY_INT; }

template<> inline int get_type_num<long>()
    { return NPY_LONG; }

template<> inline int get_type_num<float>()
    { return NPY_FLOAT; }

template<> inline int get_type_num<double>()
    { return NPY_DOUBLE; }

template<> inline int get_type_num<long double>()
    { return NPY_LONGDOUBLE; }

template<> inline int get_type_num<long long>()
    { return NPY_LONGLONG; }

/** TODO: Add the following basic Numpy types

#define NPY_SIZEOF_COMPLEX_FLOAT 8
#define NPY_SIZEOF_DOUBLE 8
#define NPY_SIZEOF_COMPLEX_DOUBLE 16
#define NPY_SIZEOF_LONGDOUBLE 16
#define NPY_SIZEOF_COMPLEX_LONGDOUBLE 32
#define NPY_SIZEOF_PY_INTPTR_T 8
#define NPY_SIZEOF_PY_LONG_LONG 8
#define NPY_SIZEOF_LONGLONG 8
*/
// --------------------------------------------------------------------

// ==============================================================
// Convert Numpy to Blitz++ Arrays

/** Makes sure a Numpy Array has the given type and dimension.  Raises
a Python exception if it does not.
@param type_num See http://docs.scipy.org/doc/numpy/reference/c-api.dtype.html eg: NPY_DOUBLE */
PyArrayObject *check_dimensions(
PyObject *ovec,
std::string const &vname,
int type_num,
int ndim,
int const *dims,
int N_ndim);


/** Convert a Numpy array to a blitz one, using the original's data (no copy).
N is the rank of the array.
@see: http://mail.scipy.org/pipermail/numpy-discussion/2004-March/002892.html

Use this template to convert from PyObject to a blitz::Array of a specific
type and dimensions.

Example:
    PyObject *Val_py;
    auto Val(py_to_blitz<double,3>(Val_py, "Val", {3, -1, 4}));

@param type_num Must match T
*/
template<class T, int N>
blitz::Array<T,N> py_to_blitz(
PyObject *ovec,
std::string const &vname,
int ndim,           // # dimensions from dimension array (user-supplied)
int const *dims)
{
    // Check data types and cast
    PyArrayObject *vec = check_dimensions(ovec, vname, get_type_num<T>(), ndim, dims, N);

    int T_size = sizeof(T);
    assert(T_size == PyArray_ITEMSIZE(vec));

    blitz::TinyVector<int,N> shape(0);
    blitz::TinyVector<int,N> strides(0);
    npy_intp *arr_dimensions = vec->dimensions;
    npy_intp *arr_strides = vec->strides;

    // Set up shape and strides
    for (int i=0;i<N;++i) {
        shape[i]   = arr_dimensions[i];
        // Python/Numpy strides are in bytes, Blitz++ in sizeof(T) units.
        strides[i] = arr_strides[i]/T_size;
    }
//printf("py_to_blitz: vec->data = %p\n", vec->data);

    return blitz::Array<T,N>((T*) PyArray_DATA(vec),shape,strides,
        blitz::neverDeleteData);
}

template<class T, int N>
blitz::Array<T,N> py_to_blitz(
PyObject *ovec,
std::string const &vname,
std::vector<int> const &dims)
    { return py_to_blitz<T,N>(ovec, vname, dims.size(), &dims[0]);
    }

// ======================================================================

template<class T>
std::vector<T> py_to_vector(PyObject *ovec, std::string const &vname, int ilen = -1)
{
    int dims[1] = {ilen};

    auto blitz_obj(py_to_blitz<T,1>(ovec, vname, 1, dims));
    int len = blitz_obj.extent(0);
    std::vector<T> ret;
    ret.reserve(len);
    for (int i=0; i<len; ++i) ret.push_back(blitz_obj(i));
    return ret;
}


template<class T>
PyObject *vector_to_py(std::vector<T> const &vec)
{
    PyObject *vec_py = NULL;
    try {
        int n = vec.size();
        int dims[] = {n};
        PyObject *vec_py = PyArray_FromDims(1, dims, get_type_num<T>());

        // Wrap as blitz for convenience
        auto vec_bl(py_to_blitz<T,1>(vec_py, "vec_py", 1, dims));

        // Copy from the vector to the Blitz array
        for (int i=0; i<n; ++i) vec_bl(i) = vec[i];
        return vec_py;
    } catch(...) {
        // De-allocated since we're not returning
        if (vec_py) Py_DECREF(vec_py);
        giss::exit(1);
    }
}

template<class T, int N>
PyObject *blitz_to_py(blitz::Array<T,N> const &array)
{
    PyObject *array_py = NULL;
    try {
        int dims[N];
        for (int i=0; i<N; ++i) dims[i] = array.extent(i);
        PyObject *array_py = PyArray_FromDims(N, dims, get_type_num<T>());

        // Wrap as blitz for convenience
        auto array_bl(py_to_blitz<T,N>(array_py, "array_py", 1, dims));

        // Copy from the arraytor to the Blitz array
        array_bl = array;
        return array_py;
    } catch(...) {
        // De-allocated since we're not returning
        if (array_py) Py_DECREF(array_py);
        giss::exit(1);
    }
}

// ======================================================================

PyObject *VectorSparseMatrix_to_py(giss::VectorSparseMatrix const &mat);

VectorSparseMatrix py_to_VectorSparseMatrix(PyObject *m_tuple, std::string const &vname);

giss::BlitzSparseMatrix py_to_BlitzSparseMatrix(PyObject *m_tuple, std::string const &vname);


}   // namespace giss
