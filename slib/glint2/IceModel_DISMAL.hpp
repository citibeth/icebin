#pragma once

#include <glint2/IceModel.hpp>
#include <glint2/Grid_XY.hpp>

namespace glint2 {

class IceModel_DISMAL : public IceModel
{
	// Assume a simple X-Y Cartesian grid
	int nx, ny;

	blitz::Array<double,2> decode(
		blitz::Array<int,1> const &indices,
		blitz::Array<double,1> const &vals);

public:

//	IceModel_DISMAL(int _nx, int _ny) : nx(_nx), ny(_ny) {}

	IceModel_DISMAL(Grid_XY const &grid) : nx(grid.nx()), ny(grid.ny()) {}

	/** Query all the ice models to figure out what fields they need */
	void get_required_fields(std::set<IceField> &fields);

	/** @param index Index of each grid value.
	@param vals The values themselves -- could be SMB, Energy, something else...
	TODO: More params need to be added.  Time, return values, etc. */
	void run_timestep(int itime,
		blitz::Array<int,1> const &indices,
		std::map<IceField, blitz::Array<double,1>> const &vals2);


};

}
