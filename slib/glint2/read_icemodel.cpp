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
	GCMCoupler const *coupler,
	NcFile &nc,
	std::string const &vname,
	IceSheet *sheet)
{
	printf("BEGIN read_icemodel(%s)\n", vname.c_str());
	auto info_var = nc.get_var((vname + ".info").c_str());
	auto const_var = nc.get_var("const");	// Physical constants

	IceModel::Type type = giss::parse_enum<IceModel::Type>(
		giss::get_att(info_var, "ice_model")->as_string(0));

	std::unique_ptr<IceModel> ice_model;
	switch(type.index()) {
		case IceModel::Type::DISMAL :
			ice_model.reset(new IceModel_DISMAL(coupler));
			break;
		case IceModel::Type::PISM :
			ice_model.reset(new glint2::gpism::IceModel_PISM(coupler));
			break;
	}

	// After this, the caller must run the following to finish IceModel setup:
	// 1. Configure the contracts
	// 2. ice_model->init(coupler, sheet->grid2, nc, vname, const_var);
	// 3. ice_model->update_ice_sheet(nc, vname, sheet);

	return ice_model;
	printf("END read_icemodel(%s)\n", vname.c_str());
}
// -----------------------------------------------------

}
