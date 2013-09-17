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

#pragma once

#include <glint2/GridDomain.hpp>
#include <unordered_set>

namespace glint2 {

// ------------------------------------------------
class GridDomain_Listed : public GridDomain {

	std::unordered_set<int> domain;
	std::unordered_set<int> halo;

public:

	void halo_add(int index_c)
	{
		halo.insert(index_c);
		domain.insert(index_c);	// Things in domain are also in halo
	}
	void domain_add(int index_c)
		{ domain.insert(index_c); }


	bool in_domain(int index_c) const
		{ return domain.find(index_c) != domain.end(); }

	bool in_halo(int index_c) const
		{ return halo.find(index_c) != halo.end(); }

	// Default implementation is OK; or re-implement to avoid
	// going through extra virtual function call
	boost::function<bool (int)> get_in_halo2();
};
// ------------------------------------------------


}
