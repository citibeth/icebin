#ifndef ICEBIN_MODELE_GRIDS_HPP
#define ICEBIN_MODELE_GRIDS_HPP

#include <icebin/modele/hntr.hpp>

namespace icebin {
namespace modele {

int const IM2 = 10800;
int const JM2 = 5400;
int const IMS = 2160;
int const JMS = 1080;
int const IMH = 720;
int const JMH = 360;
int const IM1 = 360;
int const JM1 = 180;
int const IM = 288;
int const JM = 180;

extern HntrGrid const g2mx2m;
extern HntrGrid const g10mx10m;
extern HntrGrid const ghxh;
extern HntrGrid const g1x1;
double const dLATM = 180.*60./JM;  //  latitude spacing (minutes)
extern HntrGrid const g1qx1;
extern HntrGrid const g2hx2;
extern HntrGrid const g5x4;

}}    // namespace

#endif    // guard
