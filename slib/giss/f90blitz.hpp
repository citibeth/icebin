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

#include <mpi.h>
#include <cstdio>
#include <blitz/array.h>
#include <iostream>

/** @file
C++ (Blitz) / Fortran interoperability for Assumed Shape (Fortran 90) arrays. */

namespace giss {

/**A reconstruction of the "dope vector" used to describe assumed
shape arrays in Fortran 90.  Used as parameters in C++ functions to
accept Fortran 90 arrays and re-constitute them as blitz++ arrays.
Peered with the Fortran derived types arr_spec_<x>, where <x> is a
number of dimensions.
@see f90blitz_f.py */
template<class ArrT, int rank>
struct F90Array {
	ArrT *base;				///< Beginning of the array
	ArrT *deltas[rank];		///< Used to calculate strides
	int lbounds[rank];		///< Lower bound index of each dimension (usually 1)
	int ubounds[rank];		///< Upper bound index of each dimension

	F90Array() : base(0) {}

	/** Extracts the dope vector from an existing blitz::Array.  Used
	to write C++ test code for Fortrn APIs (so we can construct and pass
	dope vectors without actually running Fortran code.) */
	F90Array(blitz::Array<ArrT,rank> &arr) {
		this->base = arr.data();

		blitz::TinyVector<int, rank> idx(0);
		for (int i=0; i<rank; ++i) {
			this->deltas[i] = this->base + arr.stride(i);
			this->lbounds[i] = arr.lbound(i);
			this->ubounds[i] = arr.ubound(i);
		}
	}

	/** Construct a blitz::Array from the dope vector.  The blitz
	array returned is identical to the original Fortran array used to
	create this dope vector. */
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
