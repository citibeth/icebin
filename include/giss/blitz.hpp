#pragma once

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

}
