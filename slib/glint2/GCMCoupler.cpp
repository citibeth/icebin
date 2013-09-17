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

#include <glint2/GCMCoupler.hpp>

namespace glint2 {

/** Query all the ice models to figure out what fields they need */
std::set<IceField> GCMCoupler::get_required_fields()
{
	std::set<IceField> ret;
	for (auto model = models.begin(); model != models.end(); ++model) {
		model->get_required_fields(ret);
	}
	return ret;
}


void GCMCoupler::read_from_netcdf(NcFile &nc, std::string const &vname,
	std::vector<std::string> const &sheet_names,
    giss::MapDict<std::string, IceSheet> const &sheets)
{
	printf("BEGIN GCMCoupler::read_from_netcdf()\n");
	int i = 0;
	for (auto name = sheet_names.begin(); name != sheet_names.end(); ++name) {
		models.insert(i, read_icemodel(nc, vname + "." + *name, sheets[*name]));
		++i;
	}
	printf("END GCMCoupler::read_from_netcdf()\n");
}


}
