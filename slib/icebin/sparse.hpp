#ifndef ICEBIN_SPARSE_HPP
#define ICEBIN_SPARSE_HPP

#include <spsparse/VectorCooArray.hpp>

namespace icebin {

	/** The sparse vector and matrix data types we'll use in IceBin. */
	typedef spsparse::VectorCooArray<long, double, 2> SparseMatrix;
	typedef spsparse::VectorCooArray<long, double, 1> SparseVector;

}
#endif
