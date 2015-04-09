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

#include <vector>
#include <blitz/array.h>

namespace glint2 {

/** Closure used to determine the elevation class of each GCM grid cell, based on its elevation.

This is built using a simple STL collection class, plus standard STL binary
search subroutines.  This class is
set up to work with an unusual array structure (see HeightMaxArray).

Typically, this class is used as a temporary decorator to a HeightMaxArray instance.

@see HeightClassifier::get_hclass()
*/
class HeightClassifier {
public :
	HeightClassifier(blitz::Array<double,1> const *_hcmax);

	/** Determines the elevation class of a grid cell.
	@param index Index of the grid cell (base=0) to determine elevation class.
	@param elevation Elevation of that grid cell. */
	int operator()(double elevation) const;

	double operator[](int ix) const { return (*hcmax)(ix); }

	int nhp() const { return hcmax->extent(0); }

private:
	blitz::Array<double,1> const *hcmax;

};	// class HeightClassifier

}	// namespace giss

