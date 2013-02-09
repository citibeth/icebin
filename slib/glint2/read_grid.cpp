#include <memory>

#include <giss/ncutil.hpp>
#include <glint2/Grid.hpp>
#include <glint2/Grid_XY.hpp>
#include <glint2/Grid_LonLat.hpp>
#include <glint2/ExchangeGrid.hpp>

namespace glint2 {

/** @param fname Name of file to load from (eg, an overlap matrix file)
@param vname Eg: "grid1" or "grid2" */
std::unique_ptr<Grid> read_grid(NcFile &nc, std::string const &vname)
{
	auto info_var = nc.get_var((vname + ".info").c_str());

	Grid::Type type = giss::parse_enum<Grid::Type>(
		giss::get_att(info_var, "type")->as_string(0));

	std::unique_ptr<Grid> grid;
	switch(type.index()) {
		case Grid::Type::XY :
			grid.reset(new Grid_XY());
			break;
		case Grid::Type::LONLAT :
			grid.reset(new Grid_LonLat());
			break;
		case Grid::Type::EXCHANGE :
			grid.reset(new ExchangeGrid());
			break;
		default :
			grid.reset(new Grid(type));
			break;
	}

	grid->read_from_netcdf(nc, vname);
	return grid;
}

}
