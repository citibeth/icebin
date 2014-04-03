#include <giss/cfnames.hpp>

namespace giss {

CFName const *get_cfname(std::string const &_id, std::string const &expected_units)
{
	CFName *cf = cf::standard_name_table()[_id];
	if ((expected_units != "") && (cf->canonical_units != expected_units)) {
		fprintf("CF Name '%s' has units '%s', but '%s' was expected.\n",
			cf->id.c_str(), cf->canonical_units.c_str(), cf->expected_units.c_str());
		throw std::exception();
	}
	return ;
}

}
