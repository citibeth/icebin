#include <giss/cfnames.hpp>
#include <giss/exit.hpp>

namespace giss {

CFName const *get_cfname(std::string const &_id, std::string const &expected_units)
{
	giss::CFName const &cfp = cf::standard_name_table().at(_id);
	if ((expected_units != "") && (cfp.canonical_units != expected_units)) {
		fprintf(stderr, "CF Name '%s' has units '%s', but '%s' was expected.\n",
			cfp.id.c_str(), cfp.canonical_units.c_str(), expected_units.c_str());
		giss::exit(1);
	}
	return &cfp;
}

}
