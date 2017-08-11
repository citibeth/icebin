#include <icebin/modele/clippers.hpp>
#include <icebin/gridgen/clippers.hpp>

namespace icebin {
namespace modele {

bool ice_sheet::clip(int zone, long index, double lon0, double lat0, double lon1, double lat1)
{
    // Is it in Greenland range?
    if (zone & GREENLAND)
        if (SphericalClip::lonlat(-74., 59., -10., 87.5,
            lon0, lat0, lon1, lat1)) return true;

    // Is it in Antarctica range?
    if (zone & ANTARCTICA)
        if (lat0 <= -60. || lat1 <= -60) return true;

    // Not in range of either ice sheet, discard
    return false;
}

}}    // namespace icebin::modele
