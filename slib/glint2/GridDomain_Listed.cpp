#include <glint2/GridDomain_Listed.hpp>

namespace glint2 {

boost::function<bool (int)> GridDomain_Listed::get_in_halo2()
	{ return boost::bind(&GridDomain_Listed::in_halo2, this, _1); }

}
