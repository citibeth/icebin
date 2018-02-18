#include <icebin/modele/grids.hpp>

using namespace std;

namespace icebin {
namespace modele {

// Grids we use (lon x lat)
HntrSpec const g1mx1m(IM1m, JM1m, 0., 1.);
HntrSpec const g2mx2m(IM2m, JM2m, 0., 2.);
HntrSpec const g10mx10m(IMS, JMS, 0., 10.);
HntrSpec const ghxh(IMH, JMH, 0., 30.);
HntrSpec const g1x1(IM1, JM1, 0., 60.);
HntrSpec const g1qx1(IM, JM, 0., dLATM);
HntrSpec const g2hx2(IMB, JMB, 0., dLATM * 2.);
HntrSpec const g5x4(72, 45, 0., dLATM * 4.);

std::map<std::string, HntrSpec const *> make_grids()
{
    std::map<std::string, HntrSpec const *> ret;

    ret.insert(make_pair(std::string("g1mx1m"), &g1mx1m));
    ret.insert(make_pair("g2mx2m", &g2mx2m));
    ret.insert(make_pair("g10mx10m", &g10mx10m));
    ret.insert(make_pair("ghxh", &ghxh));
    ret.insert(make_pair("g1x1", &g1x1));
    ret.insert(make_pair("g1qx1", &g1qx1));
    ret.insert(make_pair("g2hx2", &g2hx2));
    ret.insert(make_pair("g5x4", &g5x4));
    
    return ret; 
}

std::map<std::string, HntrSpec const *> const grids(make_grids());

}}
