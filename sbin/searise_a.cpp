#include <boost/bind.hpp>

#include <glint2/Grid_XY.hpp>
#include <glint2/clippers.hpp>
//#include <glint2/constants.hpp>
#include <boost/filesystem.hpp>

static const double km = 1000.0;

using namespace glint2;

int main(int argc, char **argv)
{
	printf("------------- Set up the local ice grid\n");

	Grid_XY grid;
	grid.name = "searise";
	grid.sproj = "+proj=stere +lon_0=0 +lat_0=-90 +lat_ts=71.0 +ellps=WGS84";

	// The true exact SeaRISE grid
	set_xy_boundaries(grid,
		(-2800.0 - 2.5)*km, (-2800.0 + 1200*5.0 + 2.5)*km,   5.0*km,
		(-2800.0 - 2.5)*km, (-2800.0 + 1200*5.0 + 2.5)*km,   5.0*km);

	grid.realize(boost::bind(&EuclidianClip::keep_all, _1));

	printf("Ice grid has %ld cells\n", grid.ncells_full());

	// ------------- Write it out to NetCDF
	fflush(stdout);
	printf("// ------------- Write it out to NetCDF\n");
	boost::filesystem::path exe_path(argv[0]);
	grid.to_netcdf(exe_path.stem().string() + ".nc");
//	grid.to_netcdf("searise.nc");
}
