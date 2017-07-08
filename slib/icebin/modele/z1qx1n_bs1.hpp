#ifndef ICEBIN_Z1QX1N_BS1_HPP
#define ICEBIN_Z1QX1N_BS1_HPP

#include <string>
#include <vector>
#include <ibmisc/blitz.hpp>
#include <ibmisc/IndexSet.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/filesystem.hpp>

//#include <cstdio>
//#include <boost/algorithm/string.hpp>
//#include <icebin/modele/hntr.hpp>


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

class HntrGrid;
extern HntrGrid const g2mx2m;
extern HntrGrid const g10mx10m;
extern HntrGrid const ghxh;
extern HntrGrid const g1x1;
double const dLATM = 180.*60./JM;  //  latitude spacing (minutes)
extern HntrGrid const g1qx1;

template<class TypeT, int RANK>
struct ArrayMeta {
    std::string const name;
    blitz::Array<TypeT, RANK> arr;
    std::array<std::string,RANK> const dims;
    std::string description;

    ArrayMeta(
        std::string const &_name,
        blitz::Array<TypeT, RANK> const &_arr,
        std::array<std::string,RANK> const &_dims)
    : name(_name), arr(_arr), dims(_dims) {}
};



/** Area of memory where a TOPO-generating procedure can place its outputs.
Should be pre-allocated before the generator is called. */
template<class TypeT, int RANK>
class ArrayBundle {
public:
    ibmisc::IndexSet<std::string> index;
    std::vector<ArrayMeta<TypeT, RANK>> data;
//    blitz::TinyVector<int,RANK> const &shape_t;

public:
//    ArrayBundle(
//        blitz::TinyVector<int,RANK> const &_shape_t) :
//        shape_t(_shape_t) {}


    std::vector<std::string> const &keys() const
        { return index.keys(); }


    /** Add a self-allocated array */
//    int add(std::string const &name);

    /** Add a self-allocated array */
    int add(
        std::string const &name,
        blitz::TinyVector<int,RANK> const &shape,
        std::array<std::string,RANK> const &dims);


    /** Add an existing array (must be Fortran-style) */
    int add(
        std::string const &name,
        blitz::Array<TypeT, RANK> &arr,
        std::array<std::string,RANK> const &dims);


    /** Add an existing array (must be Fortran-style) */
    int add(
        std::string const &name,
        blitz::Array<TypeT, RANK> arr,
        std::array<std::string,RANK> const &dims);

    ArrayMeta<TypeT, RANK> &at(std::string const &name)
        { return data[index.at(name)]; }

    void ncio(ibmisc::NcIO &ncio, std::string const &prefix, std::string const &snc_type);

};


/** Add a self-allocated array */
template<class TypeT, int RANK>
int ArrayBundle<TypeT,RANK>::add(
    std::string const &name,
    blitz::TinyVector<int, RANK> const &shape,
    std::array<std::string,RANK> const &dims)
{
    size_t ix = index.insert(name);
    // No array provided, allocate a new one
    data.push_back(ArrayMeta<TypeT,RANK>(
        name,
        blitz::Array<TypeT,RANK>(shape, blitz::fortranArray),
        dims));
    return ix;
}


/** Add an existing array (must be Fortran-style) */
template<class TypeT, int RANK>
int ArrayBundle<TypeT,RANK>::add(
    std::string const &name,
    blitz::Array<TypeT, RANK> &arr,
    std::array<std::string,RANK> const &dims)
{
    size_t ix = index.insert(name);
    // Array provided, reference it
    data.push_back(ArrayMeta<TypeT,RANK>(name, arr, dims));
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


// --------------------------------------------------------------------


template<class TypeT, int RANK>
void ArrayBundle<TypeT,RANK>::ncio(ibmisc::NcIO &ncio, std::string const &prefix, std::string const &snc_type)
{
    // NetCDF stores in row-major (C) order.  Reverse dimensions, and convert array to C-style

    for (size_t i=0; i<index.size(); ++i) {
        auto &meta(data[i]);

        auto dims_f(ibmisc::get_or_add_dims(ncio,
            meta.arr,
            ibmisc::to_vector(meta.dims)));

        std::vector<netCDF::NcDim> dims_c;
        for (int j=dims_f.size()-1; j>=0; --j) dims_c.push_back(dims_f[j]);

        auto &arr_c(ncio.tmp.copy<blitz::Array<TypeT,RANK>>(
            ibmisc::f_to_c(meta.arr)));

        auto ncvar(ibmisc::ncio_blitz(ncio, arr_c, false, prefix + meta.name, snc_type, dims_c));
        ncvar.putAtt("description", meta.description);

    }
}



/** Area of memory where a TOPO-generating procedure can place its outputs.
Should be pre-allocated before the generator is called. */
class TopoOutputs {
    ArrayBundle<double, 2> bundle;
public:
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

/** Input files:
 Z2MX2M.NGDC = FOCEN2: Ocean Fraction (0 or 1)
               ZETOP2: Solid Topography (m) except for ice shelves
    Z10MX10M = FLAKES: Lake Fraction (0:1)
     ZICEHXH = dZGICH: Glacial Ice Thickness (m)
               FGICEH: Glacial Ice Fraction (0:1)
               ZSOLDH: Ice Top Topography (m)
      ZNGDC1 = FCONT1: Continent Fraction (0:1)
               FGICE1: Glacial Ice Fraction (0:1) */
class TopoInputs {
public:
    ArrayBundle<double, 2> bundle;

    blitz::Array<double, 2> FOCEN2, ZETOP2;
    blitz::Array<double, 2> FLAKES;
    blitz::Array<double, 2> dZGICH, FGICEH, ZSOLDH;
    blitz::Array<double, 2> FCONT1, FGICE1;

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
