#ifndef ICEBIN_MODELE_CLIPPERS_HPP
#define ICEBIN_MODELE_CLIPPERS_HPP

namespace icebin {
namespace modele {

struct ice_sheet {
    static const int GREENLAND = 1;
    static const int ANTARCTICA = 2;

    /** Clipping function clips around ice sheets */
    static bool clip(int zone, double lon0, double lat0, double lon1, double lat1);
};




}}    // namespace
#endif    // guard
