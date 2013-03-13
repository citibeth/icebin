#pragma once

#include <giss/Dict.hpp>
#include <glint2/IceModel.hpp>

namespace glint2 {

class GCMCoupler {
public:

	// Only needed by root MPI node in MPI version
	giss::MapDict<int,IceModel> models;

	/** Query all the ice models to figure out what fields they need */
	std::set<IceField> get_required_fields();

	virtual void read_from_netcdf(NcFile &nc, std::string const &vname,
		std::vector<std::string> const &sheet_names);

};

}
