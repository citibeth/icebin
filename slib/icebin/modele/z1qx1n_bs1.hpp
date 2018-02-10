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


template<int RANK>
extern ibmisc::ArrayBundle<double,RANK> topo_outputs_bundle();

inline ibmisc::ArrayBundle<double,2> topo_outputs_bundle2(bool allocate)
{
    auto ret(topo_outputs_bundle<2>());
    if (allocate) ret.allocate({IM,JM}, {"im", "jm"},
            true, blitz::fortranArray);
    return ret;
}

/** Area of memory where a TOPO-generating procedure can place its outputs.
Should be pre-allocated before the generator is called. */
template <int RANK>
class TopoOutputs {
public:
    ibmisc::ArrayBundle<double,RANK> bundle;

    // 0 or 1, Bering Strait 1 cell wide         GISS 1Qx1,
    blitz::Array<double,RANK> &FOCEAN;
    // Lake Surface Fraction (0:1)                GISS 1Qx1,
    blitz::Array<double,RANK> &FLAKE;
    // Ground Surface Fraction (0:1)              GISS 1Qx1,
    blitz::Array<double,RANK> &FGRND;
    // Glacial Ice Surface Fraction except Greenland (0:1)         GISS 1Qx1,
    blitz::Array<double,RANK> &FGICE;
    // Greenland Ice surface
    blitz::Array<double,RANK> &FGICE_greenland;
    // Atmospheric Topography (m)                 ETOPO2 1Qx1,
    blitz::Array<double,RANK> &ZATMO;
    // Ocean Thickness (m)                       ETOPO2 1Qx1,
    blitz::Array<double,RANK> &dZOCEN;
    // Lake Thickness (m)                        ETOPO2 1Qx1,
    blitz::Array<double,RANK> &dZLAKE;
    // Glacial Ice Thickness (m)                 Ekholm,Bamber,
    blitz::Array<double,RANK> &dZGICE;
    // Solid Ground Topography (m)               ETOPO2 1Qx1,
    blitz::Array<double,RANK> &ZSOLDG;
    // Lowest Solid Topography (m)                ETOPO2 1Qx1,
    blitz::Array<double,RANK> &ZSGLO;
    // Lake Surface Topography (m)                ETOPO2 1Qx1,
    blitz::Array<double,RANK> &ZLAKE;
    // Topography Break between Ground and GIce   ETOPO2 1Qx1,
    blitz::Array<double,RANK> &ZGRND;
    // Highest Solid Topography (m)               ETOPO2 1Qx1/
    blitz::Array<double,RANK> &ZSGHI;
    // Fractional Ocean Cover (0:1)              ETOPO2 1Qx1/
    blitz::Array<double,RANK> &FOCENF;

    TopoOutputs(ibmisc::ArrayBundle<double,RANK> &&_bundle);

};

template<int RANK>
ibmisc::ArrayBundle<double,RANK> topo_outputs_bundle();

template<int RANK>
ibmisc::ArrayBundle<double,RANK> topo_outputs_bundle()
{
    ibmisc::ArrayBundle<double,RANK> bundle;
    bundle.add("FOCEAN", {
        "description", "0 or 1, Bering Strait 1 cell wide",
        "units", "1",
        "source", "GISS 1Qx1",
    });
    bundle.add("FLAKE", {
        "description", "Lake Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    });
    bundle.add("FGRND", {
        "description", "Ground Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    });
    bundle.add("FGICE", {
        "description", "Glacial Ice Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    });
    bundle.add("FGICE_greenland", {
        "description", "Greenland Ice Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    });
    bundle.add("ZATMO", {
        "description", "Atmospheric Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    });
    bundle.add("dZOCEN", {
        "description", "Ocean Thickness",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    });
    bundle.add("dZLAKE", {
        "description", "Lake Thickness",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    });
    bundle.add("dZGICE", {
        "description", "Glacial Ice Thickness",
        "units", "m",
        "sources", "Ekholm,Bamber",
    });
    bundle.add("ZSOLDG", {
        "description", "Solid Ground Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    });
    bundle.add("ZSGLO", {
        "description", "Lowest Solid Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    });
    bundle.add("ZLAKE", {
        "description", "Lake Surface Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    });
    bundle.add("ZGRND", {
        "description", "Topography Break between Ground and GIce",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    });
    bundle.add("ZSGHI", {
        "description", "Highest Solid Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    });
    bundle.add("FOCENF", {
        "description", "Fractional ocean ocver",
        "units", "1",
        "sources", "GISS 1Qx1",
    });
#if 0
    if (allocate) {
        bundle.allocate(blitz::shape(IM,JM), {"im", "jm"},
            true, blitz::fortranArray);
    }
#endif
    return bundle;
}

template<int RANK>
TopoOutputs<RANK>::TopoOutputs(ibmisc::ArrayBundle<double,RANK> &&_bundle) :
    bundle(std::move(_bundle)),
    FOCEAN(bundle.array("FOCEAN")),
    FLAKE(bundle.array("FLAKE")),
    FGRND(bundle.array("FGRND")),
    FGICE(bundle.array("FGICE")),
    FGICE_greenland(bundle.array("FGICE_greenland")),
    ZATMO(bundle.array("ZATMO")),
    dZOCEN(bundle.array("dZOCEN")),
    dZLAKE(bundle.array("dZLAKE")),
    dZGICE(bundle.array("dZGICE")),
    ZSOLDG(bundle.array("ZSOLDG")),
    ZSGLO(bundle.array("ZSGLO")),
    ZLAKE(bundle.array("ZLAKE")),
    ZGRND(bundle.array("ZGRND")),
    ZSGHI(bundle.array("ZSGHI")),
    FOCENF(bundle.array("FOCENF"))
{}



extern ibmisc::ArrayBundle<double,2> topo_inputs_bundle(bool allocate);

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

    TopoInputs(ibmisc::ArrayBundle<double,2> &&_bundle);
};


extern ibmisc::ArrayBundle<double,2> greenland_inputs_bundle(bool allocate);
class GreenlandInputs {
public:
    ibmisc::ArrayBundle<double, 2> bundle;

    // --- 2 minutes (IM2, JM2)
    // FOCEAN = Ocean Fraction (0:1)                     NGDC 2x2 (minutes)
    blitz::Array<double, 2> &FOCEN2;
    blitz::Array<double, 2> &ZETOP2;
    // ...regridded to 10-minute grid
    blitz::Array<double, 2> &FOCENS;

    // --- 10 minute (IMS, JMS)
    // FLAKE: Lake Fraction (0 to 1)
    blitz::Array<double, 2> &FLAKES;

    // --- 1/2-degree (IMH, JMH)
    // Ice Fraction (m) Ekholm, Bamber
    blitz::Array<double, 2> &FGICEH;

    // --- 1-degree (IM1, JM1)
    // PLAND: = 1-POCEAN
    blitz::Array<double, 2> &FCONT1;
    // PLICE: % OF LAND ICE
    blitz::Array<double, 2> &FGICE1;

    GreenlandInputs(ibmisc::ArrayBundle<double,2> &&_bundle);

};

extern void read_raw(TopoInputs &in, bool separate, GreenlandInputs *greenland, ibmisc::FileLocator const &files);

// =================================================================

/** Output from procedure to generate hi-res ice cover and elevations */
struct Etopo1Ice {
    // Outputs of Etopo1Ice() procedure
    spsparse::SparseSet<int,int> dimI;
    /** Contains: fgiceId, elevId
    	(mask is implied, since masked-out cells are not included in dimI) */
    ibmisc::ArrayBundle<double,1> bundle;

    Etopo1Ice(ibmisc::ArrayBundle<double,1> &&_bundle) :
        bundle(std::move(_bundle)) {}
};

// =================================================================

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

extern void z1qx1n_bs1(TopoInputs &in, TopoOutputs<2> &out);


}}
#endif     // guard
