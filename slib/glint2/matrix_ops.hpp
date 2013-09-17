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

#include <blitz/array.h>
#include <giss/SparseAccumulator.hpp>
#include <giss/SparseMatrix.hpp>
#include <giss/blitz.hpp>
#include <glint2/Grid.hpp>
#include <glint2/Grid_LonLat.hpp>
#include <glint2/Grid_XY.hpp>

namespace glint2 {

void divide_by(giss::VectorSparseMatrix &mat,
	giss::SparseAccumulator<int,double> &area,
	giss::SparseAccumulator<int,double> &area_inv);



} // namespace glint2

