#pragma once

namespace glint2 {

/** Utility to get height-classified indexing right.
TODO: Make this more efficient by rounding n1 up to a power of 2. */
class HCIndex {
public:

	int const n1;

	HCIndex(int _n1) : n1(_n1) {}

	int ik_to_index(int i, int k)	// k == ihc == hc
		{ return k * n1 + i; }

	int index_to_ik(int index, int &i, int &k)
	{
		k = index / n1;
		i = index - k*n1;
	}
};

}
