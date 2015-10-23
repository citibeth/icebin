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

#include <giss/Proj2.hpp>
#include <giss/ncutil.hpp>
#include <giss/exit.hpp>

namespace giss {

Proj2::Proj2(std::string const &_sproj, Direction _direction) :
sproj(_sproj), direction(_direction)
	{ realize(); }


Proj2::Proj2(Proj2 const &rhs) :
sproj(rhs.sproj), direction(rhs.direction)
	{ realize(); }


void Proj2::realize()
{
	if (sproj == "") return;

	_proj = Proj(sproj);
	_llproj = _proj.latlong_from_proj();
}

/** Transforms a single coordinate pair
@param src Source coordinate system
@param dest Destination coordinate system.
@param x0 Source x (or longitude) coordinate (radians)
@param y0 Source y (or latitude) coordinate (radians)
@param x1 Destination x (or longitude) coordinate (radians)
@param y1 Destination y (or latitude) coordinate (radians) */
int Proj2::transform(double x0, double y0, double &x1, double &y1) const
{
	if (!is_valid()) {
		x1 = x0;
		y1 = y0;
		return 0;
	}

	if (direction == Direction::XY2LL) {
		int ret = giss::transform(_proj, _llproj, x0, y0, x1, y1);
		x1 *= R2D;
		y1 *= R2D;
		return ret;
	}

	x0 *= D2R;
	y0 *= D2R;
	return giss::transform(_llproj, _proj, x0, y0, x1, y1);
}




void Proj2::netcdf_define(NcFile &nc, NcVar *info_var, std::string const &vname) const
{
	// ------ Attributes
	if (is_valid()) {
		info_var->add_att(vname.c_str(), sproj.c_str());
		std::string sdir = (direction == Direction::XY2LL ? "xy2ll" : "ll2xy");
		info_var->add_att((vname + ".direction").c_str(), sdir.c_str());
	} else {
		info_var->add_att(vname.c_str(), "");
		info_var->add_att((vname + ".direction").c_str(), "");
	}
}


/** @param fname Name of file to load from (eg, an overlap matrix file)
@param vname Eg: "grid1" or "grid2" */
void Proj2::read_from_netcdf(
NcFile &nc,
NcVar *info_var,
std::string const &vname)
{
	sproj = std::string(get_att(info_var, vname.c_str())->as_string(0));
	if (sproj == "") {
		clear();
	} else {
		std::string sdir(std::string(get_att(info_var, (vname + ".direction").c_str())->as_string(0)));
		if (sdir == "xy2ll") direction = Direction::XY2LL;
		else if (sdir == "ll2xy") direction = Direction::LL2XY;
		else giss::exit(1);
	}
}


}	// namespace giss
