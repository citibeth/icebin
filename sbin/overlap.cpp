#include <boost/bind.hpp>

#include <glint2/Grid.hpp>
#include <boost/filesystem.hpp>
#include <netcdfcpp.h>
#include <string>
#include <glint2/ExchangeGrid.hpp>

static const double km = 1000.0;

using namespace glint2;

int main(int argc, char **argv)
{
	std::string fname1(argv[1]);
	std::string fname2(argv[2]);

	printf("------------- Set up the projection\n");
	double proj_lon_0 = -39;
	double proj_lat_0 = 90;
	char sproj[100];
	sprintf(sproj,
		"+proj=stere +lon_0=%f +lat_0=%f +lat_ts=71.0 +ellps=WGS84",
		proj_lon_0, proj_lat_0);


	printf("------------- Read grid1\n");
	NcFile nc1(fname1.c_str());
	auto grid1(glint2::read_grid(nc1, "grid"));
	nc1.close();

	printf("------------- Read grid2\n");
	NcFile nc2(fname2.c_str());
	auto grid2(glint2::read_grid(nc2, "grid"));
	nc2.close();

	printf("--------------- Overlapping\n");
	ExchangeGrid exch(*grid1, *grid2);

	printf("--------------- Writing Out\n");
	std::string fname = grid1->name + "-" + grid2->name + ".nc";
//	NcFile nc(fname.c_str(), NcFile::Replace);
	exch.to_netcdf(fname);
}
