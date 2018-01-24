#include <icebin/modele/grids.hpp>

namespace icebin {
namespace modele {

// Grids we use (lon x lat)
HntrSpec const g2mx2m(IM2, JM2, 0., 2.);
HntrSpec const g10mx10m(IMS, JMS, 0., 10.);
HntrSpec const ghxh(IMH, JMH, 0., 30.);
HntrSpec const g1x1(IM1, JM1, 0., 60.);
HntrSpec const g1qx1(IM, JM, 0., dLATM);
HntrSpec const g2hx2(144, 90, 0., dLATM * 2.);
HntrSpec const g5x4(72, 45, 0., dLATM * 4.);

}}
