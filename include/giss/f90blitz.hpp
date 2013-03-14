#include <cstdio>
#include <blitz/array.h>

namespace giss {

/** Used to accept Fortran 90 arrays and re-constitute them as blitz++ arrays.
@see f90blitz_f.py */
template<class ArrT, int rank>
struct F90Array {
	ArrT *base;
	ArrT *deltas[rank];
	int lbounds[rank];
	int ubounds[rank];

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

}	// namespace giss
