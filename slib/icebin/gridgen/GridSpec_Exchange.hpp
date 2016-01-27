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

#include <memory>
#include <icebin/Grid.hpp>
#include <ibmisc/Proj.hpp>

namespace icebin {

class GridSpec_Exchange : public GridSpec {
public :
	Grid const *gridA;
	Grid const *gridI;

	/** Transformation to get from local Grid coords to Exchange Grid coords */
	std::string sproj;

	/** Keep track of the "full" indexing space for the Overlap Matrix. */
	long gridA_cells_nfull;
	long gridI_cells_nfull;

	GridSpec_Exchange() : gridA(0), gridI(0), gridA_cells_nfull(-1), gridI_cells_nfull(-1) {}

	void make_grid(Grid &grid);
	void ncio(ibmisc::NcIO &ncio, std::string const &vname);
};



}	// namespace glint2
