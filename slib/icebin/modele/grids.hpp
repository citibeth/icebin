#ifndef ICEBIN_MODELE_GRIDS_HPP
#define ICEBIN_MODELE_GRIDS_HPP

#include <icebin/modele/hntr.hpp>

namespace icebin {
namespace modele {

int const IM1m = 21600;
int const JM1m = 10801;
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

extern HntrSpec const g1mx1m;
extern HntrSpec const g2mx2m;
extern HntrSpec const g10mx10m;
extern HntrSpec const ghxh;
extern HntrSpec const g1x1;
double const dLATM = 180.*60./JM;  //  latitude spacing (minutes)
extern HntrSpec const g1qx1;
extern HntrSpec const g2hx2;
extern HntrSpec const g5x4;

}}    // namespace

#endif    // guard
