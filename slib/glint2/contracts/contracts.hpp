#pragma once

#include <string>

namespace glint2 {

namespace contracts {

	/** Indicates the grid this field is supposed to be on. */
	const unsigned GRID_BITS = 3;

	const unsigned ATMOSPHERE = 1;
	const unsigned ICE = 2;
	const unsigned ELEVATION = 3;

	/** Field is returned at initialization time, before the first coupling. */
	const unsigned INITIAL = 4;

	extern std::string to_str(unsigned int flags);

};



}
