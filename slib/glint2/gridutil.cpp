#include <snowdrift/maputils.hpp>
#include <snowdrift/geometry.hpp>
#include <snowdrift/Proj.hpp>
#include <snowdrift/constants.hpp>

namespace glint2 {

// ---------------------------------------------------------
/** Project a latitude line segment
n >=1, number of points to add to polygoon.  Don't add (lat,lon1). */
void LineProjector::proj_latitude(Cell &cell, int n,
	double lon0, double lon1, double lat)
{

	for (int i=0; i<n; ++i) {
		double lon = lon0 + (lon1-lon0) * ((double)i/(double)n);
		double x,y;
		int err = transform(llproj, proj, lon*D2R, lat*D2R, x, y);
		add_vertex(cell, x, y);
	}
}

/** Project a latitude line segment
n >=1, number of points to add to polygoon.  Don't add (lat,lon1). */
void LineProjector::proj_meridian(Cell &cell, int n,
	double lon, double lat0, double lat1)
{
	for (int i=0; i<n; ++i) {
		double lat = lat0 + (lat1-lat0) * ((double)i/(double)n);
		double x,y;
		transform(llproj, proj, lon*D2R, lat*D2R, x, y);
		// proj.ll2xy(lon, lat, x, y);
		add_vertex(cell, x, y);
	}
}

}
