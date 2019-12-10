#ifndef ICEBIN_MODELE_GRIDS_HPP
#define ICEBIN_MODELE_GRIDS_HPP

#include <icebin/GridSpec.hpp>

namespace icebin {
namespace modele {

int const IM1m = 21600;
int const JM1m = 10800;
int const JM5m = 4320;
int const IM5m = 2160;
int const JM1m_gridreg = 10801;    // Grid-registered ETOPO1
int const IM2m = 10800;
int const JM2m = 5400;
int const IMS = 2160;        // 10-minute
int const JMS = 1080;
int const IMH = 720;
int const JMH = 360;
int const IM1 = 360;
int const JM1 = 180;
int const IM = 288;
int const JM = 180;
int const IMB = 144;     // g2hx2
int const JMB = 90;


extern HntrSpec const g1mx1m;
extern HntrSpec const g2mx2m;
extern HntrSpec const g10mx10m;
extern HntrSpec const gqxq;
extern HntrSpec const ghxh;
extern HntrSpec const g1x1;
double const dLATM = 180.*60./JM;  //  latitude spacing (minutes)
extern HntrSpec const g1qx1;
extern HntrSpec const g2hx2;
extern HntrSpec const g5x4;

extern std::map<std::string, HntrSpec const *> const grids;

const double EQ_RAD = 6.371e6; /// Radius of the Earth (same as in ModelE)

// See LIGrid.F90
const int UI_UNUSED = 0;
const int UI_ICEBIN = 1;
const int UI_NOTHING= 2;


}}    // namespace icebin::modele

#endif    // guard
