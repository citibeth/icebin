#include <giss/Proj2.hpp>
#include <giss/ncutil.hpp>

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
		else throw std::exception();
	}
}


}	// namespace giss
