#pragma once

#include <mpi.h>
#include <cstdio>
#include <blitz/array.h>
#include <iostream>

namespace giss {

/** Used to accept Fortran 90 arrays and re-constitute them as blitz++ arrays.
@see f90blitz_f.py */
template<class ArrT, int rank>
struct F90Array {
	ArrT *base;
	ArrT *deltas[rank];
	int lbounds[rank];
	int ubounds[rank];

	/** Extract F90Array info from existing blitz::Array.  Used to
	write C++ test code for Fortrn APIs. */
	F90Array(blitz::Array<ArrT,rank> &arr) {
		this->base = arr.data();

		blitz::TinyVector<int, rank> idx(0);
		for (int i=0; i<rank; ++i) {
			this->deltas[i] = this->base + arr.stride(i);
			this->lbounds[i] = arr.lbound(i);
			this->ubounds[i] = arr.ubound(i);
		}
	}


	blitz::Array<ArrT,rank> to_blitz()
	{
		blitz::TinyVector<int, rank> shape, stride;
		blitz::GeneralArrayStorage<rank> stor;
		for (int i=0; i<rank; ++i) {
			shape[i] = ubounds[i] - lbounds[i] + 1;
			stride[i] = deltas[i] - base;
			stor.base()[i] = lbounds[i];
			// Ordering is not needed because we're using stride
			// stor.ordering()[i] = i;      // Fortran ordering, blitz++ manual p31
		}
		return blitz::Array<ArrT,rank>(base, shape, stride,
			blitz::neverDeleteData, stor);
	}

};

/** Print details of array */
template<class ArrT, int rank>
std::ostream& operator<<(std::ostream& os, F90Array<ArrT, rank> const &arr_f)
{
	char buf[32];
	sprintf(buf, "%p", arr_f.base);

	os << "F90Array(" << std::string(buf) << " [";
	for (int i=0; i<rank; ++i) os << (arr_f.deltas[i] - arr_f.base) << " ";
	os << "] : [";
	for (int i=0; i<rank; ++i) os << arr_f.lbounds[i] << " ";
	os << "] -- [";
	for (int i=0; i<rank; ++i) os << arr_f.ubounds[i] << " ";
	os << "])";
    return os;
}


}	// namespace giss
