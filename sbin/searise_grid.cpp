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
	grid.sproj = "+proj=stere +lon_0=-39 +lat_0=90 +lat_ts=71.0 +ellps=WGS84";

#define CHOICE 3

#if CHOICE==1
	set_xy_boundaries(grid,
		0., 10., 5.,
		0., 10., 5.);
#endif

#if CHOICE==2
	// The true exact SeaRISE grid
	set_xy_boundaries(grid,
		(- 800.0 - 2.5)*km, (- 800.0 + 300.0*5 + 2.5)*km,   5*km,
		(-3400.0 - 2.5)*km, (-3400.0 + 560.0*5 + 2.5)*km,   5*km);
#endif

#if CHOICE==3
	// Approximate 50km SeaRISE-like grid
	set_xy_boundaries(grid,
		(-800)*km, (-800 + 300*5)*km,     50*km,
		(-3400)*km, (-3400 + 560*5)*km,   50*km);
#endif

	grid.realize(boost::bind(&EuclidianClip::keep_all, _1));

	printf("Ice grid has %ld cells\n", grid.ncells_full);

	// ------------- Write it out to NetCDF
	fflush(stdout);
	printf("// ------------- Write it out to NetCDF\n");
	boost::filesystem::path exe_path(argv[0]);
//	grid.to_netcdf(exe_path.stem().string() + ".nc");
	grid.to_netcdf("searise.nc");
}
