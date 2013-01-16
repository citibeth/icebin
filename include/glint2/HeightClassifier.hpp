#pragma once

#include <vector>
#include <blitz/array.h>

namespace giss {

/** Closure used to determine the elevation class of each GCM grid cell, based on its elevation.

This is built using a simple STL collection class, plus standard STL binary
search subroutines.  This class is
set up to work with an unusual array structure (see HeightMaxArray).

Typically, this class is used as a temporary decorator to a HeightMaxArray instance.

@see HeightClassifier::get_hclass()
*/
class HeightClassifier {
public :
	HeightClassifier(blitz::Array<double,1> *_hcmax) : hcmax(_hcmax) {}

	/** Determines the elevation class of a grid cell.
	@param index Index of the grid cell (base=0) to determine elevation class.
	@param elevation Elevation of that grid cell. */
	int operator()(double elevation);

	double operator[](int ix) { return (*hcmax)(ix); }

private:
	blitz::Array<double,1> *hcmax;

};	// class HeightClassifier

}	// namespace giss

