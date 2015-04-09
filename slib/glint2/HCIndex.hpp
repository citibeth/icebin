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
		(MODELE)	(1)		// ModelE-style indexing: C++ (nhp, n1)
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
