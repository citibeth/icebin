#include <cmath>
#include <ibmisc/blitz.hpp>

namespace icebin {
namespace modele {

struct HntrGrid {
    int im;    // Number of cells in east-west direction
    int jm;    // Number of cells in north-south direction
    blitz::Array<double,1> dxyp;

    int size() { return im * jm; }

    // number (fraction) of cells in east-west direction from
    // International Date Line (180) to western edge of cell IA=1
    double offi;

    // minutes of latitude for non-polar cells on grid A
    double dlat;

    HntrGrid(int _im, int _jm, double _offi, double _dlat) :
        im(_im), jm(_jm), offi(_offi), dlat(_dlat), dxyp(jm)
    {
        // Calculate the sperical area of grid cells
        double dLON = 2.*M_PI / im;
        double dLAT = M_PI / jm;
        for (int j=1; j<=jm; ++j) {
            double SINS = sin(dLAT*(J-jm/2-1));
            double SINN = sin(dLAT*(J-jm/2));
            dxyp(J) = dLON * (SINN - SINS)
        }
    }

    template<class TypeT>
    blitz::Array<TypeT, 2> Array()
        { return blitz::Array(im,jm, blitz::fortranArray); }
};

class Hntr {
    HntrGrid const Agrid;
    HntrGrid const Bgrid;

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
    blitz::Array<int,1> IMIN, IMAX
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

    Input: WTA = weighting array for values on the A grid
             A = per unit area or per unit mass quantity
    Output:  B = horizontally interpolated quantity on B grid
    */
    void regrid(
        blitz::Array<double,1> const &WTA,
        blitz::Array<double,1> const &A,
        blitz::Array<double,1> &B);

    /** HNTR4P is similar to HNTR4 but polar values are replaced by their
    longitudinal mean. */
    void regridp(
        blitz::Array<double,1> const &WTA,
        blitz::Array<double,1> const &A,
        blitz::Array<double,1> &B);

private:
    void partition_east_west();
    void partition_north_south();

};    // class Hntr

}}

