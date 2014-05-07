#include <giss/CouplingContract.hpp>
#include <ostream>

namespace giss {

int CouplingContract::add_field(CoupledField &&cf)
{
	int index = _ix_to_field.size();
	_name_to_ix.insert(std::make_pair(cf.name, index));
	bool is_unit = (cf.name == "unit");
	_ix_to_field.push_back(std::move(cf));
	if (is_unit) {
		_unit_ix = index;
	} else {
		_size_nounit += 1;
	}
	return index;
}


std::ostream &operator<<(std::ostream &out, CouplingContract const &con) {
	for (auto field = con.begin(); field != con.end(); ++field) {
		out << "    " << *field << std::endl;
	}
	return out;
}

}
