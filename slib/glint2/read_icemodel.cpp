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

#include <glint2/pism/IceModel_PISM.hpp>	// PISM insists on being included first
#include <glint2/IceModel.hpp>
#include <glint2/IceModel_DISMAL.hpp>

namespace glint2 {

/** @param sheet (OPTIONAL): Info on the ice sheet data structure */
std::unique_ptr<IceModel> read_icemodel(
	MPI_Comm gcm_comm,
	boost::filesystem::path const &config_dir,		/** Directory of glint2 configuration file */
	NcFile &nc,
	std::string const &vname,
	IceSheet const *sheet)
{
	printf("BEGIN read_icemodel(%s)\n", vname.c_str());
	auto info_var = nc.get_var((vname + ".info").c_str());
	auto const_var = nc.get_var("const");	// Physical constants

	IceModel::Type type = giss::parse_enum<IceModel::Type>(
		giss::get_att(info_var, "ice_model")->as_string(0));

	std::unique_ptr<IceModel> ice_model;
	switch(type.index()) {
		case IceModel::Type::DISMAL : {
			Grid_XY const *grid2 = dynamic_cast<Grid_XY const *>(&*(sheet->grid2));
			auto dismal_var = nc.get_var((vname + ".dismal").c_str());	// DISMAL parameters
			ice_model.reset(new IceModel_DISMAL(*grid2, config_dir, dismal_var, const_var));
		} break;
		case IceModel::Type::PISM : {
			std::shared_ptr<Grid_XY const> grid2 = std::dynamic_pointer_cast<Grid_XY const>(sheet->grid2);
//			Grid_XY const *grid2 = dynamic_cast<Grid_XY const *>(&*(sheet->grid2));
			auto pism_var = nc.get_var((vname + ".pism").c_str());	// PISM parameters
			ice_model.reset(new glint2::pism::IceModel_PISM(grid2, gcm_comm, config_dir, pism_var, const_var));
		} break;
#if 0
		case IceModel::Type::ISSM :
			ice_model.reset(new IceModel_ISSM(const_var));
			break;
#endif
	}

	ice_model->read_from_netcdf(nc, vname);
	return ice_model;
	printf("END read_icemodel(%s)\n", vname.c_str());
}
// -----------------------------------------------------

}
