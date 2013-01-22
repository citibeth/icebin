#include <boost/bind.hpp>

#include <glint2/Grid_XY.hpp>
#include <glint2/clippers.hpp>
//#include <glint2/constants.hpp>
#include <boost/filesystem.hpp>

static const double km = 1000.0;

using namespace glint2;

int main(int argc, char **argv)
{
	// ------------- Set up the local ice grid
	printf("// ------------- Set up the local ice grid\n");

	Grid_XY grid;
	grid.name = "searise";


	set_xy_boundaries(grid,
		0., 10., 5.,
		0., 10., 5.);

#if 0
#if 1
	// The true exact SeaRISE grid
	set_xy_boundaries(grid,
		(- 800.0 - 2.5)*km, (- 800.0 + 300.0*5 + 2.5)*km,   5*km,
		(-3400.0 - 2.5)*km, (-3400.0 + 560.0*5 + 2.5)*km,   5*km);
#else
	// Approximate 50km SeaRISE-like grid
	set_xy_boundaries(grid,
		(-800)*km, (-800 + 300*5)*km,     50*km,
		(-3400)*km, (-3400 + 560*5)*km,   50*km);
#endif
#endif

	grid.realize(boost::bind(&EuclidianClip::keep_all, _1));

	printf("Ice grid has %ld cells\n", grid.ncells_full);

	// ------------- Write it out to NetCDF
	fflush(stdout);
	printf("// ------------- Write it out to NetCDF\n");
	boost::filesystem::path exe_path(argv[0]);
	grid.to_netcdf(exe_path.stem().string() + ".nc");
}
