#pragma once

#include <boost/enum.hpp>

namespace giss {

template<class T>
inline T parse_enum(char const *str) {
	auto ret = T::get_by_name(str);
	if (!ret) {
		fprintf(stderr, "Error converting from string '%s' for boost::enum type %s\n", str, typeid(T).name());
		throw std::exception();
	}
	return *ret;
}

}
