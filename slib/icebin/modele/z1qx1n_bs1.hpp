#ifndef ICEBIN_Z1QX1N_BS1_HPP
#define ICEBIN_Z1QX1N_BS1_HPP

#include <string>
#include <vector>
#include <ibmisc/blitz.hpp>
#include <ibmisc/IndexSet.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/filesystem.hpp>
#include <icebin/modele/grids.hpp>

namespace icebin {
namespace modele {


/** Array meta-data stored with GISS-type files */
template<class TypeT, int RANK>
struct ArrayMeta {
    std::string const name;
    blitz::Array<TypeT, RANK> arr;
    std::array<std::string,RANK> const dims;
    std::string description;
    std::string units;
    std::string source;

    ArrayMeta(
        std::string const &_name,
        blitz::Array<TypeT, RANK> const &_arr,
        std::array<std::string,RANK> const &_dims,
        std::string const &_description = "",
        std::string const &_units = "",
        std::string const &_source = "")

    : name(_name), arr(_arr), dims(_dims), description(_description), units(_units), source(_source) {}
};



/** Area of memory where a TOPO-generating procedure can place its outputs.
Should be pre-allocated before the generator is called. */
template<class TypeT, int RANK>
class ArrayBundle {
public:
    ibmisc::IndexSet<std::string> index;
    std::vector<ArrayMeta<TypeT, RANK>> data;

public:
    ArrayMeta<TypeT, RANK> const &at(std::string const &name) const
        { return data[index.at(name)]; }


    std::vector<std::string> const &keys() const
        { return index.keys(); }

    /** Add a self-allocated array */
    int add(
        std::string const &name,
        blitz::TinyVector<int,RANK> const &shape,
        std::array<std::string,RANK> const &dims,
        std::string const &description = "",
        std::string const &units = "",
        std::string const &source = "");

    /** Add an existing array (must be Fortran-style) */
    int add(
        std::string const &name,
        blitz::Array<TypeT, RANK> arr,
        std::array<std::string,RANK> const &dims);

    int add(ArrayMeta<TypeT, RANK> const &datum)
    {
        size_t ix = index.insert(datum.name);
        data.push_back(datum);
    }

    ArrayMeta<TypeT, RANK> &at(std::string const &name)
        { return data[index.at(name)]; }

    void ncio(ibmisc::NcIO &ncio, std::string const &prefix, std::string const &snc_type, std::vector<std::string> const &vars = {"<all>"});

};


/** Add a self-allocated array */
template<class TypeT, int RANK>
int ArrayBundle<TypeT,RANK>::add(
    std::string const &name,
    blitz::TinyVector<int, RANK> const &shape,
    std::array<std::string,RANK> const &dims,
    std::string const &description,
    std::string const &units,
    std::string const &source)
{
    size_t ix = index.insert(name);
    // No array provided, allocate a new one
    data.push_back(ArrayMeta<TypeT,RANK>(
        name,
        blitz::Array<TypeT,RANK>(shape, blitz::fortranArray),
        dims, description, units, source));
    return ix;
}

/** Add an existing array (must be Fortran-style) */
template<class TypeT, int RANK>
int ArrayBundle<TypeT,RANK>::add(
    std::string const &name,
    blitz::Array<TypeT, RANK> arr,
    std::array<std::string,RANK> const &dims)
{
    size_t ix = index.insert(name);
    // Array provided, reference it
    data.push_back(ArrayMeta<TypeT,RANK>(name, arr, dims));
    return ix;
}

// --------------------------------------------------------------------

template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::ncio(ibmisc::NcIO &ncio, std::string const &prefix, std::string const &snc_type, std::vector<std::string> const &vars)
{
    std::vector<std::string> all_vars;
    std::vector<std::string> const *myvars;
    if (vars.size() == 1 && vars[0] == "<all>") {
        for (size_t i=0; i<index.size(); ++i) {
            all_vars.push_back(index[i]);
        }
        myvars = &all_vars;
    } else {
        myvars = &vars;
    }

    for (auto &var : *myvars) {
        int i=index.at(var);

        auto &meta(data[i]);

        // Set up the dimensions (Fortran order)
        auto dims_f(ibmisc::get_or_add_dims(ncio,
            meta.arr,
            ibmisc::to_vector(meta.dims)));

        // Write the NetCDF variable
        // (will auto-reverse dims if it detects column major)
        auto ncvar(ibmisc::ncio_blitz(ncio, meta.arr, false, prefix + meta.name, snc_type, dims_f));
        ncvar.putAtt("description", meta.description);
        ncvar.putAtt("units", meta.units);
        ncvar.putAtt("source", meta.source);

    }
}



/** Area of memory where a TOPO-generating procedure can place its outputs.
Should be pre-allocated before the generator is called. */
class TopoOutputs {
public:
    ArrayBundle<double, 2> bundle;

    // 0 or 1, Bering Strait 1 cell wide         GISS 1Qx1,
    blitz::Array<double, 2> FOCEAN;
    // Lake Surface Fraction (0:1)                GISS 1Qx1,
    blitz::Array<double, 2> FLAKE;
    // Ground Surface Fraction (0:1)              GISS 1Qx1,
    blitz::Array<double, 2> FGRND;
    // Glacial Ice Surface Fraction (0:1)         GISS 1Qx1,
    blitz::Array<double, 2> FGICE;
    // Atmospheric Topography (m)                 ETOPO2 1Qx1,
    blitz::Array<double, 2> ZATMO;
    // Ocean Thickness (m)                       ETOPO2 1Qx1,
    blitz::Array<double, 2> dZOCEN;
    // Lake Thickness (m)                        ETOPO2 1Qx1,
    blitz::Array<double, 2> dZLAKE;
    // Glacial Ice Thickness (m)                 Ekholm,Bamber,
    blitz::Array<double, 2> dZGICE;
    // Solid Ground Topography (m)               ETOPO2 1Qx1,
    blitz::Array<double, 2> ZSOLDG;
    // Lowest Solid Topography (m)                ETOPO2 1Qx1,
    blitz::Array<double, 2> ZSGLO;
    // Lake Surface Topography (m)                ETOPO2 1Qx1,
    blitz::Array<double, 2> ZLAKE;
    // Topography Break between Ground and GIce   ETOPO2 1Qx1,
    blitz::Array<double, 2> ZGRND;
    // Highest Solid Topography (m)               ETOPO2 1Qx1/
    blitz::Array<double, 2> ZSGHI;
    // Fractional Ocean Cover (0:1)              ETOPO2 1Qx1/
    blitz::Array<double, 2> FOCENF;

    TopoOutputs(ArrayBundle<double, 2> &&_bundle);
};

TopoOutputs make_topo_outputs();

/** Input files:
 Z2MX2M.NGDC = FOCEN2: Ocean Fraction (0 or 1)
               ZETOP2: Solid Topography (m) except for ice shelves
    Z10MX10M = FLAKES: Lake Fraction (0:1)

Ice Sheets:
     ZICEHXH = dZGICH: Glacial Ice Thickness (m)
               FGICEH: Glacial Ice Fraction (0:1)
               ZSOLDH: Ice Top Topography (m)

Mountain Glaciers:
      ZNGDC1 = FCONT1: Continent Fraction (0:1); ALL continents, incl. Greenland & Antarctica
               FGICE1: Glacial Ice Fraction (0:1)
*/
class TopoInputs {
public:
    ArrayBundle<double, 2> bundle;

    // --- 2 minutes (IM2, JM2)
    // FOCEAN = Ocean Fraction (0:1)                     NGDC 2x2 (minutes)
    blitz::Array<double, 2> FOCEN2;
    // ZSOLID = Solid Ground Topography (m)              NGDC 2x2 (minutes)
    blitz::Array<double, 2> ZETOP2;

    // --- 10 minute (IMS, JMS)
    // FLAKE: Lake Fraction (0 to 1)
    blitz::Array<double, 2> FLAKES;

    // --- 1/2-degree (IMH, JMH)
    // Ice Thickness (m) Ekholm, Bamber
    blitz::Array<double, 2> dZGICH;
    // Ice Fraction (m) Ekholm, Bamber
    blitz::Array<double, 2> FGICEH;
    // Ice Top Topography (m) ETOP05, Bamber
    blitz::Array<double, 2> ZSOLDH;

    // --- 1-degree (IM1, JM1)
    // PLAND: = 1-POCEAN
    blitz::Array<double, 2> FCONT1;
//    // Greenland separated out from original FCONT1
//    blitz::Array<double, 2> fcont1_greenland;


    // PLICE: % OF LAND ICE
    blitz::Array<double, 2> FGICE1;

    TopoInputs(ArrayBundle<double, 2> &&_bundle);
};

/** Allocates our own TopoInputs */
TopoInputs make_topo_inputs();


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
