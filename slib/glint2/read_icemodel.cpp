/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <glint2/IceModel.hpp>
#include <glint2/IceModel_DISMAL.hpp>

namespace glint2 {

/** @param sheet (OPTIONAL): Info on the ice sheet data structure */
std::unique_ptr<IceModel> read_icemodel(NcFile &nc, std::string const &vname, IceSheet const *sheet)
{
	printf("BEGIN read_icemodel(%s)\n", vname.c_str());
	auto info_var = nc.get_var((vname + ".info").c_str());

	IceModel::Type type = giss::parse_enum<IceModel::Type>(
		giss::get_att(info_var, "ice_model")->as_string(0));

	std::unique_ptr<IceModel> ice_model;
	switch(type.index()) {
		case IceModel::Type::DISMAL : {
			Grid_XY const *grid2 = dynamic_cast<Grid_XY const *>(&*(sheet->grid2));
			ice_model.reset(new IceModel_DISMAL(*grid2));
		} break;
#if 0
		case IceModel::Type::PISM :
			ice_model.reset(new IceModel_PISM());
			break;
		case IceModel::Type::ISSM :
			ice_model.reset(new IceModel_ISSM());
			break;
#endif
	}

	ice_model->read_from_netcdf(nc, vname);
	return ice_model;
	printf("END read_icemodel(%s)\n", vname.c_str());
}
// -----------------------------------------------------

}
