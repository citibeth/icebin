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

