//#include <snowdrift/maputils.hpp>
#include <glint2/gridutil.hpp>
#include <giss/Proj.hpp>
#include <giss/constant.hpp>
#include <boost/bind.hpp>

namespace glint2 {


Vertex *VertexCache::add_vertex(double x, double y) {
	auto ii = _vertices.find(std::make_pair(x,y));
	if (ii != _vertices.end()) {	// Found in cache
		return ii->second;
	} else {
		return grid->add_vertex(Vertex(x,y));
	}
}

Vertex *VertexCache::add_vertex(Cell &cell, double x, double y)
{
	Vertex *v = add_vertex(x,y);	// Add to the Grid
	cell.add_vertex(v);
	return v;
}

// ======================================================
// LineProjector

LineProjector::LineProjector(giss::Proj const &_proj, giss::Proj const &_llproj,
	boost::function<Vertex * (double, double)> const &add_vertex) :
	proj(_proj), llproj(_llproj), _add_vertex(add_vertex)
{}

/** Convenience constructor. */
LineProjector::LineProjector(giss::Proj const &_proj, VertexCache *vcache) :
	proj(_proj),
	llproj(proj.latlong_from_proj()),
	_add_vertex(boost::bind(&VertexCache::add_vertex, vcache, _1, _2))
{}

// ---------------------------------------------------------
/** Project a latitude line segment
n >=1, number of points to add to polygoon.  Don't add (lat,lon1). */
void LineProjector::proj_latitude(Cell &cell, int n,
	double lon0, double lon1, double lat)
{

	for (int i=0; i<n; ++i) {
		double lon = lon0 + (lon1-lon0) * ((double)i/(double)n);
		double x,y;
		int err = transform(llproj, proj, lon*giss::D2R, lat*giss::D2R, x, y);
		cell.add_vertex(_add_vertex(x, y));
	}
}

/** Project a latitude line segment
n >=1, number of points to add to polygoon.  Don't add (lat,lon1). */
void LineProjector::proj_meridian(Cell &cell, int n,
	double lon, double lat0, double lat1)
{
	for (int i=0; i<n; ++i) {
		double n_inv = 1.0 / (double)n;
		double lat = lat0 + (lat1-lat0) * (double)i * n_inv;
		double x,y;
		transform(llproj, proj, lon*giss::D2R, lat*giss::D2R, x, y);
		// proj.ll2xy(lon, lat, x, y);
		cell.add_vertex(_add_vertex(x, y));
	}
}

}
