#pragma once

#include <glint2/Grid_LonLat.hpp>

namespace glint2 {
namespace giss {

/** @param lons As read out of ModelE netCDF file */
extern void set_lonlat_centers(glint2::Grid_LatLon &grid,
	std::vector<double> const &lons,
	std::vector<double> const &lats);

extern void set_lonlat_4x5(glint2::Grid_LatLon &grid);
extern void set_lonlat_2x2_5(glint2::Grid_LatLon &grid);


}}
