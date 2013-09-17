/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <glint2/HeightClassifier.hpp>
#include <algorithm>

namespace glint2 {

HeightClassifier::HeightClassifier(
	blitz::Array<double,1> const *_hcmax) :
	hcmax(_hcmax) {}

int HeightClassifier::operator()(double elevation) const
{
//	end = hcmax.extent(0) - 1;	// Last height class always ends at infinity

	auto begin(hcmax->begin());
	auto end(hcmax->end());
	--end;		// Last height class always ends at infinity
	auto top(std::upper_bound(begin, end, elevation));

	int hc = top.position()[0] - begin.position()[0];
//printf("hc = %d\n", hc);
	return hc;
}

}
