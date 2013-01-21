#include <glint2/Grid_LonLat.hpp>

namespace glint2 {
namespace giss {

/** @param lons As read out of ModelE netCDF file */
void set_lonlat_centers_giss(glint2::Grid_LonLat &grid,
	std::vector<double> const &lons,
	std::vector<double> const &lats)
{
	// --------- Reprocess lat/lon format for a quadrilateral mesh
	// (Assume latlon grid)
	// Shift lats to represent edges of grid boxes

	// Set up latitudes, splitting the difference between grid cell centers
	// and taking care of polar caps
	{size_t n = lats.size();
		grid.latb.clear();
		grid.latb.reserve(n-1);
		grid.latb.push_back(lats[1] - (lats[2] - lats[1])*.5);
		for (int j=1; j < n-2; ++j) {
			grid.latb.push_back((lats[j] + lats[j+1]) * .5);
		}
		grid.latb.push_back(lats[n-2] + (lats[n-1] - lats[n-2]) * .5);
		grid.south_pole = true;
		grid.north_pole = true;
	}

	// Shift lons to represent edges of grid boxes
	{size_t n = lons.size();
		grid.lonb.clear();
		grid.lonb.reserve(n+1);
		grid.lonb.push_back(.5 * (lons[0] + lons[n-1] - 360.));	// Assume no overlap
		for (int i=1; i < n; ++i) {
			grid.lonb.push_back(.5 * (lons[i] + lons[i-1]));
		}
		// Do repeat first point at the end.
		grid.lonb.push_back(grid.lonb[0] + 360.);
	}

	grid.south_pole = true;
	grid.north_pole = true;
}

// -------------------------------------------------------------------
// Latitude and longitude gridcell boundaries for the 4x5 grid

static const std::vector<double> lonb_4x5 = {-180,-175,-170,-165,-160,-155,-150,-145,-140,-135,-130,-125,-120,-115,-110,-105,-100,-95,-90,-85,-80,-75,-70,-65,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180};

static const std::vector<double> latb_4x5 = {-88,-84,-80,-76,-72,-68,-64,-60,-56,-52,-48,-44,-40,-36,-32,-28,-24,-20,-16,-12,-8,-4,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88};

/** Just initializes lonb and latb, doesn't do full set-up */
void set_lonlat_4x5(glint2::Grid_LonLat &grid)
{
	grid.lonb = lonb_4x5;
	grid.latb = latb_4x5;
	grid.south_pole = true;
	grid.north_pole = true;
}

void set_lonlat_2x2_5(glint2::Grid_LonLat &grid)
{
	// Create the 2x2.5 grid from the 4x5 grid.
	grid.lonb.clear();
	for (int i=0; i<lonb_4x5.size()-1; ++i) {
		grid.lonb.push_back(lonb_4x5[i]);
		grid.lonb.push_back(.5*(lonb_4x5[i] + lonb_4x5[i+1]));
	}
	grid.lonb.push_back(lonb_4x5.back());

	grid.latb.clear();
	for (int i=0; i<latb_4x5.size()-1; ++i) {
		grid.latb.push_back(latb_4x5[i]);
		grid.latb.push_back(.5*(latb_4x5[i] + latb_4x5[i+1]));
	}
	grid.latb.push_back(latb_4x5.back());

	grid.south_pole = true;
	grid.north_pole = true;
}


}}
