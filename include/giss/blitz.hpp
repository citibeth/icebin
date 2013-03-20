#pragma once

#include <cstdio>
#include <vector>
#include <blitz/array.h>

namespace giss {

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

}
