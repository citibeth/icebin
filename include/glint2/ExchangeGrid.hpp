#pragma once

#include <memory>
#include <glint2/Grid.hpp>
#include <giss/Proj.hpp>

namespace glint2 {

class ExchangeGrid : public Grid {
public :
	/** Transformation to get from local Grid coords to Exchange Grid coords */
	giss::Proj2 proj1, proj2;

	/** @param proj Projection to use to project Lon/Lat grids to XY,
	if no projection is found in the XY-type grid. */
	ExchangeGrid(Grid const &grid1, Grid const &grid2, std::string const &_sproj="");

	virtual boost::function<void()> netcdf_define(NcFile &nc, std::string const &vname) const;

	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);

};


/** @param grid2 Put in an RTree */
extern std::unique_ptr<Grid> compute_exchange_grid(Grid &grid1, Grid &grid2);

}	// namespace glint2
