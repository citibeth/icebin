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

#pragma once

#include <giss/Proj.hpp>
#include <giss/constant.hpp>
#include <netcdfcpp.h>
#include <giss/Proj.hpp>

namespace giss {


class Proj2 {
public:
	std::string sproj;
	enum class Direction {LL2XY, XY2LL};
	Direction direction;
protected:
	Proj _proj, _llproj;
	void realize();
public:

	bool is_valid() const { return _proj.is_valid(); }

	Proj2(std::string const &_sproj, Direction _direction);

	Proj2() : direction(Direction::LL2XY) {}

	void clear() {
		_proj.clear();
		_llproj.clear();
	}

	Proj2(Proj2 const &rhs);

	Proj2(Proj2 const &rhs, Direction _direction) :
		sproj(rhs.sproj), direction(_direction)
		{ realize(); }

	void init(std::string const &_sproj, Direction _direction)
	{
		sproj = _sproj;
		direction = _direction;
		realize();
	}


	/** Transforms a single coordinate pair
	@param src Source coordinate system
	@param dest Destination coordinate system.
	@param x0 Source x (or longitude) coordinate (radians)
	@param y0 Source y (or latitude) coordinate (radians)
	@param x1 Destination x (or longitude) coordinate (radians)
	@param y1 Destination y (or latitude) coordinate (radians) */
	int transform(double x0, double y0, double &x1, double &y1) const;

	void netcdf_define(NcFile &nc, NcVar *info_var, std::string const &vname) const;
	void read_from_netcdf(NcFile &nc, NcVar *info_var, std::string const &vname);
	
};

}
