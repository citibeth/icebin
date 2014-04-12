#include <glint2/CouplingContract.hpp>
#include <ostream>

namespace glint2 {

std::ostream &operator<<(std::ostream &out, CouplingContract const &con) {
	for (auto ii = con._ix_to_field.begin(); ii != con._ix_to_field.end(); ++ii) {
		CoupledField const &cpf(*ii);
		out << "    ";
		out << cpf;
		out << std::endl;
	}
	return out;
}

}
