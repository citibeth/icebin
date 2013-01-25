#pragma once

#include <cassert>
#include <vector>
#include <exception>
#include <blitz/array.h>
#include <giss/SparseMatrix.hpp>

namespace giss {

PyObject *init_module(
std::string const &module_name,
std::string const &module_description,
PyMethodDef **function_sets,	// NULL-terminated Array of NULL-terminated arrays
PyTypeObject **types);			// NULL-terminated Array of type pointers

// ====================================================================
// --------------------------------------------------------------------
// Convert template types to Numpy type_nums

template<class T>
inline int get_type_num()
{
	PyErr_SetString(PyExc_ValueError, "get_type_num(): Unknown type_num");
	throw std::exception();
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
This should be called AFTER appropriate typechecking has taken place.
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
int ndim,			// # dimensions from dimension array (user-supplied)
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
    return blitz::Array<T,N>((T*) vec->data,shape,strides,
		blitz::neverDeleteData);
}

template<class T, int N>
blitz::Array<T,N> py_to_blitz(
PyObject *ovec,
std::string const &vname,
std::vector<int> const &dims)
	{ return py_to_blitz<T,N>(ovec, vname, dims.size(),
		(int)dims.size(), &dims[0]);
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
		throw std::exception();
	}
}

// ======================================================================

PyObject *VectorSparseMatrix_to_py(giss::VectorSparseMatrix const &mat);

VectorSparseMatrix py_to_VectorSparseMatrix(PyObject *m_tuple, std::string const &vname);

}	// namespace giss
