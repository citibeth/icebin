#include <giss/IndexTranslator.hpp>

namespace giss {

/** @param n Size of space a (runs [0...n-1])
@param used Set of indices that are used in space a */
void IndexTranslator::init(int size_a, std::set<int> const &used)
{
printf("IndexTranslator::init(%s, size_a=%d, size_b=%ld)\n", _name.c_str(), size_a, used.size());
	_a2b.clear(); _a2b.resize(size_a, -1);
	_b2a.clear(); _b2a.reserve(used.size());
	for (auto ia_ptr = used.begin(); ia_ptr != used.end(); ++ia_ptr) {
		int ib = _b2a.size();
		_b2a.push_back(*ia_ptr);
		_a2b[*ia_ptr] = ib;
	}
}


}
