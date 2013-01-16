#include <snowdrift/HeightClassifier.hpp>
#include <algorithm>

namespace giss {

HeightClassifier::HeightClassifier(
	std::vector<blitz::Array<double,1>> *_height_max) :
height_max(_height_max) {}

int HeightClassifier::operator()(double elevation)
{
	auto begin(hcmax.begin());
	auto end(hcmax.end());
	end -= 1;		// Last height class always ends at infinity
	iterator top(std::upper_bound(begin, end, elevation));

	int hc = top - begin;
	return hc;
}

}
