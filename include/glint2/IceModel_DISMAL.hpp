#pragma once

#include <glint2/IceModel.hpp>

namespace glint2 {

class IceModel_DISMAL : public IceModel
{

	/** Query all the ice models to figure out what fields they need */
	void get_required_fields(std::set<IceField> &fields);

	/** @param index Index of each grid value.
	@param vals The values themselves -- could be SMB, Energy, something else...
	TODO: More params need to be added.  Time, return values, etc. */
	void run_timestep(
		blitz::Array<int,1> const &indices,
		std::map<IceField, blitz::Array<double,1>> const &vals2);


};

}
