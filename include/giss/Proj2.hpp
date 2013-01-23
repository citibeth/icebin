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


	int ll2xy(double lon0, double lat0, double &x1, double &y1) const
	{
		lon0 *= D2R;
		lat0 *= D2R;
		return giss::transform(_llproj, _proj, lon0, lat0, x1, y1);
	}

	int xy2ll(double x0, double y0, double &lon1, double &lat1) const
	{
		int ret = giss::transform(_proj, _llproj, x0, y0, lon1, lat1);
		lon1 *= R2D;
		lat1 *= R2D;
		return ret;
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
