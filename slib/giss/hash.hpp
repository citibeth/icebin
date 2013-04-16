#pragma once

#include <functional>

namespace giss {

template<class A, class B>
struct HashPair {
	size_t operator()(std::pair<A, B> const &x) const throw()
		{ return std::hash<A>()(x.first)*31 + std::hash<B>()(x.second); }
};

}
