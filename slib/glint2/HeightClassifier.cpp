#include <glint2/HeightClassifier.hpp>
#include <algorithm>

namespace glint2 {

HeightClassifier::HeightClassifier(
	blitz::Array<double,1> *_hcmax) :
	hcmax(_hcmax) {}

int HeightClassifier::operator()(double elevation)
{
//	end = hcmax.extent(0) - 1;	// Last height class always ends at infinity

	auto begin(hcmax->begin());
	auto end(hcmax->end());
	--end;		// Last height class always ends at infinity
	auto top(std::upper_bound(begin, end, elevation));

	int hc = top.position()[0] - begin.position()[0];
	return hc;
}

}
