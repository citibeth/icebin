#include <memory>

#include <glint2/Grid.hpp>
#include <glint2/Grid_XY.hpp>
#include <glint2/Grid_LonLat.hpp>

namespace glint2 {

/** @param fname Name of file to load from (eg, an overlap matrix file)
@param vname Eg: "grid1" or "grid2" */
std::unique_ptr<Grid> read_grid(NcFile &nc, std::string const &vname)
{
	auto info_var = nc.get_var((vname + ".info").c_str());
	std::string stype(info_var->get_att("type")->as_string(0));

	std::unique_ptr<Grid> grid;
	if (stype == "xy") {
		grid.reset(new Grid_XY());
	} else if (stype == "lonlat") {
		grid.reset(new Grid_LonLat());
	} else {
		grid.reset(new Grid("generic"));
	}

	grid->read_from_netcdf(nc, vname);
	return grid;
}

}
