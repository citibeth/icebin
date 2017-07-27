#ifndef ICEBIN_Z1QX1N_BS1_HPP
#define ICEBIN_Z1QX1N_BS1_HPP

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <ibmisc/blitz.hpp>
#include <ibmisc/IndexSet.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/filesystem.hpp>
#include <ibmisc/bundle.hpp>
#include <icebin/modele/grids.hpp>

namespace icebin {
namespace modele {


/** Area of memory where a TOPO-generating procedure can place its outputs.
Should be pre-allocated before the generator is called. */
class TopoOutputs {
public:
    ibmisc::ArrayBundle<double,2> bundle;

    // 0 or 1, Bering Strait 1 cell wide         GISS 1Qx1,
    blitz::Array<double, 2> &FOCEAN;
    // Lake Surface Fraction (0:1)                GISS 1Qx1,
    blitz::Array<double, 2> &FLAKE;
    // Ground Surface Fraction (0:1)              GISS 1Qx1,
    blitz::Array<double, 2> &FGRND;
    // Glacial Ice Surface Fraction (0:1)         GISS 1Qx1,
    blitz::Array<double, 2> &FGICE;
    // Atmospheric Topography (m)                 ETOPO2 1Qx1,
    blitz::Array<double, 2> &ZATMO;
    // Ocean Thickness (m)                       ETOPO2 1Qx1,
    blitz::Array<double, 2> &dZOCEN;
    // Lake Thickness (m)                        ETOPO2 1Qx1,
    blitz::Array<double, 2> &dZLAKE;
    // Glacial Ice Thickness (m)                 Ekholm,Bamber,
    blitz::Array<double, 2> &dZGICE;
    // Solid Ground Topography (m)               ETOPO2 1Qx1,
    blitz::Array<double, 2> &ZSOLDG;
    // Lowest Solid Topography (m)                ETOPO2 1Qx1,
    blitz::Array<double, 2> &ZSGLO;
    // Lake Surface Topography (m)                ETOPO2 1Qx1,
    blitz::Array<double, 2> &ZLAKE;
    // Topography Break between Ground and GIce   ETOPO2 1Qx1,
    blitz::Array<double, 2> &ZGRND;
    // Highest Solid Topography (m)               ETOPO2 1Qx1/
    blitz::Array<double, 2> &ZSGHI;
    // Fractional Ocean Cover (0:1)              ETOPO2 1Qx1/
    blitz::Array<double, 2> &FOCENF;

    TopoOutputs(bool allocate=true);
};

/** Input files:
 Z2MX2M.NGDC = FOCEN2: Ocean Fraction (0 or 1)
               ZETOP2: Solid Topography (m) except for ice shelves
    Z10MX10M = FLAKES: Lake Fraction (0:1)

Ice Sheets:
     ZICEHXH = dZGICH: Glacial Ice Thickness (m)
               FGICEH: Glacial Ice Fraction (0:1)
               ZSOLDH: Ice Top Topography (m)

Mountain Glaciers:
      ZNGDC1 = FCONT1: Continent Fraction (0:1)
               FGICE1: Glacial Ice Fraction (0:1)
*/
class TopoInputs {
public:
    ibmisc::ArrayBundle<double, 2> bundle;

    // --- 2 minutes (IM2, JM2)
    // FOCEAN = Ocean Fraction (0:1)                     NGDC 2x2 (minutes)
    blitz::Array<double, 2> &FOCEN2;
    // ZSOLID = Solid Ground Topography (m)              NGDC 2x2 (minutes)
    blitz::Array<double, 2> &ZETOP2;

    // --- 10 minute (IMS, JMS)
    // FLAKE: Lake Fraction (0 to 1)
    blitz::Array<double, 2> &FLAKES;

    // --- 1/2-degree (IMH, JMH)
    // Ice Thickness (m) Ekholm, Bamber
    blitz::Array<double, 2> &dZGICH;
    // Ice Fraction (m) Ekholm, Bamber
    blitz::Array<double, 2> &FGICEH;
    // Ice Top Topography (m) ETOP05, Bamber
    blitz::Array<double, 2> &ZSOLDH;

    // --- 1-degree (IM1, JM1)
    // PLAND: = 1-POCEAN
    blitz::Array<double, 2> &FCONT1;
    // PLICE: % OF LAND ICE
    blitz::Array<double, 2> &FGICE1;

    TopoInputs(bool allocate=true);
};

// =================================================================

extern void read_raw(TopoInputs &in, ibmisc::FileLocator const &files);


/*
Z1QX1N.BS1.F    Create Z1QX1N, Nonfractional ocean    2011/11/15

Fortran Compile and Go:  FCG Z1QX1N.BS1.CRE HNTR4.o

Z1QX1N.BS1.CRE creates a DataFile of observed surface fractions,
topography, and water and ice thicknesses.

Check the following:
  "1QX1"
  "IM  ="
  "JM  ="
  "dLATM =" assumes equal spacing in latitude including poles
  "HNTR40"
  "FILOUT"
  Lines following "Set following grid cells to be continent"
  Lines following "Set following grid cells to be ocean"

Input files:
 Z2MX2M.NGDC = FOCEN2: Ocean Fraction (0 or 1)
               ZETOP2: Solid Topography (m) except for ice shelves
    Z10MX10M = FLAKES: Lake Fraction (0:1)
     ZICEHXH = dZGICH: Glacial Ice Thickness (m)
               FGICEH: Glacial Ice Fraction (0:1)
               ZSOLDH: Ice Top Topography (m)
      ZNGDC1 = FCONT1: Continent Fraction (0:1)
               FGICE1: Glacial Ice Fraction (0:1)
*/
extern void callZ(
    // (IM2, JM2)
    blitz::Array<double,2> &FOCEN2,
    blitz::Array<double,2> &ZSOLD2,
    blitz::Array<double,2> &ZSOLG2,

    // (IM, IM)
    blitz::Array<double,2> &FOCEAN,
    blitz::Array<double,2> &FLAKE,
    blitz::Array<double,2> &FGRND,

    // (IM, IM)
    blitz::Array<double,2> &ZATMO,
    blitz::Array<double,2> &dZLAKE,
    blitz::Array<double,2> &ZSOLDG,
    blitz::Array<double,2> &ZSGLO,
    blitz::Array<double,2> &ZLAKE,
    blitz::Array<double,2> &ZGRND,
    blitz::Array<double,2> &ZSGHI);

extern void z1qx1n_bs1(TopoInputs &in, TopoOutputs &out);


}}
#endif     // guard
