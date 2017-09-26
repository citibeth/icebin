#include <icebin/modele/grids.hpp>

namespace icebin {
namespace modele {

// Grids we use (lon x lat)
HntrGrid const g2mx2m(IM2, JM2, 0., 2.);
HntrGrid const g10mx10m(IMS, JMS, 0., 10.);
HntrGrid const ghxh(IMH, JMH, 0., 30.);
HntrGrid const g1x1(IM1, JM1, 0., 60.);
HntrGrid const g1qx1(IM, JM, 0., dLATM);
HntrGrid const g2hx2(144, 90, 0., dLATM * 2.);
HntrGrid const g5x4(72, 45, 0., dLATM * 4.);

}}
