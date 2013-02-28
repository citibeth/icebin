#include <glint2/GridDomain.hpp>

namespace glint2 {

boost::function<bool (int)> GridDomain::get_in_halo()
{ return boost::bind(this, &GridDomain::in_halo, _1); }

}
