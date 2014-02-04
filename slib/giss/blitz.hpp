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

#include <cstdio>
#include <vector>
#include <blitz/array.h>

/** \file blitz.hpp
Type conversion utilities related to Blitz++ arrays.  Blitz++ arrays
are so flexible... if subroutines are written to take them as
arguments, then other array data types can be converted as needed. */

namespace giss {

/** Converts a std::vector to a Blitz++ 1-D array that shares the same memory.
@param vec The input vector */
template<class T>
blitz::Array<T,1> vector_to_blitz(std::vector<T> &vec)
{
    blitz::TinyVector<int,1> shape(0);
    blitz::TinyVector<int,1> strides(0);

	shape[0] = vec.size();
	strides[0] = 1;		// Blitz++ strides in sizeof(T) units

    return blitz::Array<T,1>(&vec[0], shape, strides,
		blitz::neverDeleteData);
}

/** Converts a const std::vector to a const Blitz++ 1-D array that shares the same memory. */
template<class T>
blitz::Array<T,1> const vector_to_blitz(std::vector<T> const &vec)
{
    blitz::TinyVector<int,1> shape(0);
    blitz::TinyVector<int,1> strides(0);

	shape[0] = vec.size();
	strides[0] = 1;		// Blitz++ strides in sizeof(T) units

	// const_cast because Blitz++ can't construct a const Array
	T *vecp = const_cast<T *>(&vec[0]);
    return blitz::Array<T,1>(vecp, shape, strides,
		blitz::neverDeleteData);
}

/** Makes sure a blitz::Array dimension.
Raises a Python exception if it does not. */

/** Checks that the dimensions of an array are what we think they
should be.  This is used in Python interface code to add sophisticated
type checking to functions.

@param vname Name of this variables (used in error comments)
@param arr The array to check dimensions
@param dims The expected dimensions.  If a dims[i] < 0, then dimension i is not checked. */
template<class T, int rank>
void check_dimensions(
std::string const &vname,
blitz::Array<T, rank> const &arr,
std::vector<int> const &dims)
{
	for (int i=0; i<rank; ++i) {
		if (dims[i] >= 0 && arr.extent(i) != dims[i]) {
			fprintf(stderr,
				"Error in %s: expected dimension #%d = %d (is %d instead)\n",
				vname.c_str(), i, dims[i], arr.extent(i));
			throw std::exception();
		}
	}
}
// ------------------------------------------------------------

/** Changes a C-style Blitz++ array (biggest stride in first dimension
and zero-based indexing) to a Fortran-style Blitz++ array (biggest
stride in last dimension and one-based indexing) that shares the same memory. */
template<class T, int rank>
blitz::Array<T, rank> c_to_f(blitz::Array<T, rank> &arr);

template<class T, int rank>
blitz::Array<T, rank> c_to_f(blitz::Array<T, rank> &arr)
{
	// Initialize an 11-dim vector of transpositions
	// (because transpose() doesn't take a TinyVector)
	int const max_dims = 11;
	int rev[max_dims];
	for (int i=rank; i<max_dims; ++i) rev[i] = 0;

	// Reverse dimensions
	for (int i=0; i<rank; ++i) rev[i] = rank-i-1;
	auto ret(arr.transpose(rev[0], rev[1], rev[2], rev[3], rev[4], rev[5], rev[6], rev[7], rev[8], rev[9], rev[10]));

	// Re-base to 1
	blitz::TinyVector<int, rank> base(1);
	ret.reindexSelf(base);

	return ret;
}
// ------------------------------------------------------------
/** Changes a C-style Blitz++ array (biggest stride in first dimension
and zero-based indexing) to a Fortran-style Blitz++ array (biggest
stride in last dimension and one-based indexing) that shares the same memory. */
template<class T, int rank>
blitz::Array<T, rank> f_to_c(blitz::Array<T, rank> &arr);

template<class T, int rank>
blitz::Array<T, rank> f_to_c(blitz::Array<T, rank> &arr)
{
	// Initialize an 11-dim vector of transpositions
	// (because transpose() doesn't take a TinyVector)
	int const max_dims = 11;
	int rev[max_dims];
	for (int i=rank; i<max_dims; ++i) rev[i] = 0;

	// Reverse dimensions
	for (int i=0; i<rank; ++i) rev[i] = rank-1-i;
	auto ret(arr.transpose(rev[0], rev[1], rev[2], rev[3], rev[4], rev[5], rev[6], rev[7], rev[8], rev[9], rev[10]));

	// Re-base to 0
	blitz::TinyVector<int, rank> base(0);
	ret.reindexSelf(base);

	return ret;
}


#define RESHAPE_BODY \
	/* Check dimensions */ \
	long src_n = 1; \
	for (int i=0; i<src_ndim; ++i) src_n *= src.extent[i]; \
	long dest_n = 1; \
	for (int i=0; i<dest_ndim; ++i) dest_n *= dest_shape[i]; \
	if (src_n != dest_n) { \
		fprintf(stderr, "blitz.hpp, giss::reshape(): Total dimension mismatch, src=%ld, dest=%ld\n", src_n, dest_n); \
		throw std::exception(); \
	} \
 \
	/* Do the reshaping */ \
	return blitz::Array<T,dest_ndim>(src.data(), dest_shape, blitz::neverDeleteData)



/** Reshape an array.  As long as src and dest have same total number
of elements.  Assumes a dense array on both sides. */
template<class T, int src_ndim, int dest_ndim>
extern blitz::Array<T, dest_ndim> reshape(
	blitz::Array<T, src_ndim> &src,
	blitz::TinyVector<int,dest_ndim> const dest_shape)
{ RESHAPE_BODY; }

/** Reshape an array.  As long as src and dest have same total number
of elements.  Assumes a dense array on both sides. */
template<class T, int src_ndim, int dest_ndim>
extern blitz::Array<T, dest_ndim> const reshape(
	blitz::Array<T, src_ndim> const &src,
	blitz::TinyVector<int,dest_ndim> const dest_shape)
{ RESHAPE_BODY; }

}
