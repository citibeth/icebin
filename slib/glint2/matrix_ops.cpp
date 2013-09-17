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

#include <cstdio>
#include <giss/Proj.hpp>
#include <glint2/matrix_ops.hpp>
#include <glint2/util.hpp>
#include <glint2/HCIndex.hpp>
#include <giss/constant.hpp>

namespace glint2 {

/** Divides mat /= area.
@param area (IN) Vector to divide by.  Value is moved out of this.
@param area_inv (OUT) 1/area */
void divide_by(giss::VectorSparseMatrix &mat,
	giss::SparseAccumulator<int,double> &area,
	giss::SparseAccumulator<int,double> &area_inv)
{
	// Compute 1 / area
	for (auto ii = area.begin(); ii != area.end(); ++ii)
		ii->second = 1.0d / ii->second;
	area_inv = std::move(area);

	// Divide by area.  (Now area is really equal to 1/area)
	for (auto ii = mat.begin(); ii != mat.end(); ++ii)
		ii.val() *= area_inv[ii.row()];
}




} // namespace glint2
