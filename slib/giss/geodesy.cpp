#include <giss/constant.hpp>

namespace giss {

inline double sqr(double x) { return x*x; }

/** See: http://www.cs.nyu.edu/visual/home/proj/tiger/gisfaq.html
@return Distance (in degrees) */
extern double haversine_distance(
double lon1_deg, double lat1_deg,
double lon2_deg, double lat2_deg)
{
        // Convert inputs to degrees
        double lon1 = lon1_deg * D2R;
        double lat1 = lat1_deg * D2R;
        double lon2 = lon2_deg * D2R;
        double lat2 = lat2_deg * D2R;

        // Apply the Haversine Formula
        double dlon = lon2 - lon1;
        double dlat = lat2 - lat1;
        double a = sqr(sin(dlat/2)) + cos(lat1) * cos(lat2) * sqr(sin(dlon/2));
        double c = 2 * atan2( sqrt(a), sqrt(1-a) );
        return c * R2D;         // Convert to degrees
}

}
