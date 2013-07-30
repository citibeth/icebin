#pragma once

#include <boost/enum.hpp>
#include <memory>
#include <cstring>

namespace glint2 {

class MatrixMaker;

/** Utility to get height-classified indexing right.
TODO: Make this more efficient by rounding n1 up to a power of 2. */
class HCIndex {
public:

	BOOST_ENUM_VALUES( Type, int,
		(UNKNOWN)	(0)
		(MODELE)	(1)		// ModelE-style indexing: C++ (nhc, n1)
	)

	virtual ~HCIndex() {}

	/** @param i Horizontal (X-Y) combined index (i1).
	To be broken down further, depending on atmosphere indexing scheme.
	@param k Vertical (elevation point) index */
	virtual int ik_to_index(int i, int k) const = 0;

	/** @param i3 Index in elevation grid */
	virtual void index_to_ik(int i3, int &i, int &k) const = 0;

	/** Factory method to instantiate a new HCIndex */
	static std::unique_ptr<HCIndex> new_HCIndex(
		Type const type,
		MatrixMaker const &mm);
};


}
