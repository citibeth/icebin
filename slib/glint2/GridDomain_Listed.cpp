#include <glint2/GridDomain.hpp>

namespace glint2 {

boost::function<bool (int)> GridDomain_Listed::get_in_halo()
{ return boost::bind(this, &GridDomain_Listed::in_halo, _1); }

}
