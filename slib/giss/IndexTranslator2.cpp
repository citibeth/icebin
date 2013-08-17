#include <giss/IndexTranslator2.hpp>

namespace giss {

/** @param n Size of space a (runs [0...n-1])
@param used_a Set of indices that are used in space a */
void IndexTranslator2::init(std::map<int, size_t> *size_a, std::set<std::pair<int,int>> const &used_a)
{

printf("IndexTranslator2::init(%s, size_a=%d, size_b=%ld)\n", _name.c_str(), size_a->begin()->second, used_a.size());

	_size_a = size_a;
	_a2b.clear();
	_b2a.clear(); _b2a.reserve(used_a.size());
	for (auto ia_ptr = used_a.begin(); ia_ptr != used_a.end(); ++ia_ptr) {
		int ib = _b2a.size();
		_b2a.push_back(*ia_ptr);
		_a2b[*ia_ptr] = ib;
	}
}


int IndexTranslator2::a2b(std::pair<int,int> a, bool check_result) const {
	int aindex = a.first;
	if (aindex < 0 || aindex >= nindex()) {
		fprintf(stderr, "%s: aindex=%d is out of range (%d, %d)\n", _name.c_str(), aindex, 0, nindex());
		throw std::exception();
	}

	int agrid = a.second;
	if (agrid < 0 || agrid >= na(aindex)) {
		fprintf(stderr, "%s: a=%d is out of range (%d, %d)\n", _name.c_str(), a, 0, na(aindex));
		throw std::exception();
	}
	auto ib = _a2b.find(a);
	if (check_result && ib == _a2b.end()) {
		fprintf(stderr, "%s: a=(%d,%d) not found\n", _name.c_str(), a.first, a.second);
		throw std::exception();
	}
	return ib->second;
}


std::pair<int,int> IndexTranslator2::b2a(int b, bool check_result) const {
	if (b < 0 || b >= nb()) {
		fprintf(stderr, "%s: b=%d is out of range (%d, %d)\n", _name.c_str(), b, 0, nb());
		throw std::exception();
	}
	std::pair<int,int> a = _b2a[b];
	if (check_result && a.first < 0) {
		fprintf(stderr, "%s: b=%d produces invalid a=(%d,%d)\n", _name.c_str(), b, a.first,a.second);
		throw std::exception();
	}
	return a;
}



}
