#include <glint2/contracts/contracts.hpp>

namespace glint2 {

	std::string contracts::to_str(unsigned int flags)
	{
		std::string ret = "";
		switch(flags & GRID_BITS) {
			case ATMOSPHERE :
				ret += "ATMOSPHERE|";
				break;
			case ICE:
				ret += "ICE|";
				break;
			case ELEVATION:
				ret += "ELEVATION|";
				break;
			default: ;
		}

		if (flags & INITIAL) ret += "INITIAL|";

		if (ret.size() > 0) ret.resize(ret.size()-1);
		return ret;
	}

}
