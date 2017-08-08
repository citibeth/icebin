#ifndef ICEBIN_HNTR_HPP
#define ICEBIN_HNTR_HPP

#include <ibmisc/blitz.hpp>
#include <spsparse/eigen.hpp>

namespace icebin {
namespace modele {

class HntrGrid {
public:
    int const im;    // Number of cells in east-west direction
    int const jm;    // Number of cells in north-south direction

    // number (fraction) of cells in east-west direction from
    // International Date Line (180) to western edge of cell IA=1
    double const offi;

    // minutes of latitude for non-polar cells on grid A
    double const dlat;

//protected:
    blitz::Array<double,1> _dxyp;

public:
    int size() const { return im * jm; }
    double dxyp(int j) const { return _dxyp(j); }

    HntrGrid(int _im, int _jm, double _offi, double _dlat);

    template<class TypeT>
    blitz::Array<TypeT, 2> Array() const
        { return blitz::Array<TypeT,2>(im,jm, blitz::fortranArray); }

};

/** Pre-computed overlap details needed to regrid from one lat/lon
    grid to another on the sphere. */
class Hntr {
public:
    HntrGrid const Agrid;
    HntrGrid const Bgrid;

//protected:
    // SINA(JA) = sine of latitude of northern edge of cell JA on grid A
    // SINB(JB) = sine of latitude of northern edge of cell JB on grid B
    // FMIN(IB) = fraction of cell IMIN(IB) on grid A west of cell IB
    // FMAX(IB) = fraction of cell IMAX(IB) on grid A east of cell IB
    // GMIN(JB) = fraction of cell JMIN(JB) on grid A south of cell JB
    // GMAX(JB) = fraction of cell JMAX(JB) on grid A north of cell JB
    // IMIN(IB) = western most cell of grid A that intersects cell IB
    // IMAX(IB) = eastern most cell of grid A that intersects cell IB
    // JMIN(JB) = southern most cell of grid A that intersects cell JB
    // JMAX(JB) = northern most cell of grid A that intersects cell JB

    blitz::Array<double, 1> SINA, SINB;
    blitz::Array<double, 1> FMIN, FMAX;
    blitz::Array<int,1> IMIN, IMAX;
    blitz::Array<double, 1> GMIN, GMAX;
    blitz::Array<int,1> JMIN, JMAX;

    // DATMIS = missing data value inserted in output array B when
    // cell (IB,JB) has integrated value 0 of WTA
    double DATMIS;

public:
    /** Initialize overlap data structures, get ready to re-grid. */
    Hntr(HntrGrid const &_A, HntrGrid const &_B, double _DATMIS);

    /**
    HNTR4 performs a horizontal interpolation of per unit area or per
    unit mass quantities defined on grid A, calculating the quantity
    on grid B.  B grid values that cannot be calculated because the
    covering A grid boxes have WTA = 0, are set to the value DATMIS.
    The area weighted integral of the quantity is conserved.
    The 3 Real input values are expected to be Real*4.

    ** NOTE **
        All arrays use 1-based (Fortran-style) indexing!!!

    Inputs to this method must all be 1-D 1-based (Fortran-style) arrays.
    See regrid() for a method accepting "natural" 2-D 1-based arrays.

    Input: WTA = weighting array for values on the A grid
             A = per unit area or per unit mass quantity
             mean_polar: polar values are replaced by their
                 longitudinal mean.
    Output:  B = horizontally interpolated quantity on B grid
    */
    void regrid1(
        blitz::Array<double,1> const &WTA,
        blitz::Array<double,1> const &A,
        blitz::Array<double,1> &B,
        bool mean_polar = false) const;

    void matrix(
        spsparse::TupleList<int,double,2> &accum,        // The output (sparse) matrix; 0-based sparse indexing
        blitz::Array<double,1> const &_WTA);


    /** Works with 0-based or 1-based N-dimensional arrays */
    template<int RANK>
    void regrid(
        blitz::Array<double,RANK> const &WTA,
        blitz::Array<double,RANK> const &A,
        blitz::Array<double,RANK> &B,
        bool mean_polar = false) const;

    template<int RANK>
    blitz::Array<double,RANK> regrid(
        blitz::Array<double,RANK> const &WTA,
        blitz::Array<double,RANK> const &A,
        bool mean_polar = false) const;

private:
    void partition_east_west();
    void partition_north_south();

};    // class Hntr



/** Works with 0-based or 1-based N-dimensional arrays */
template<int RANK>
void Hntr::regrid(
    blitz::Array<double,RANK> const &WTA,
    blitz::Array<double,RANK> const &A,
    blitz::Array<double,RANK> &B,
    bool mean_polar) const
{
    auto B1(ibmisc::reshape1(B, 1));
    return regrid1(
        ibmisc::reshape1(WTA, 1),
        ibmisc::reshape1(A, 1),
        B1, mean_polar);
}

template<int RANK>
blitz::Array<double,RANK> Hntr::regrid(
    blitz::Array<double,RANK> const &WTA,
    blitz::Array<double,RANK> const &A,
    bool mean_polar) const
{
    blitz::Array<double,2> B(Bgrid.Array<double>());
    regrid(WTA, A, B, mean_polar);
    return B;
}


}}

#endif    // guard
