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
