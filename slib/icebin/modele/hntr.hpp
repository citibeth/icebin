#ifndef ICEBIN_HNTR_HPP
#define ICEBIN_HNTR_HPP

#include <ibmisc/blitz.hpp>
#include <icebin/eigen_types.hpp>

namespace icebin {
namespace modele {

class HntrGrid {
public:
    int im;    // Number of cells in east-west direction
    int jm;    // Number of cells in north-south direction

    // number (fraction) of cells in east-west direction from
    // International Date Line (180) to western edge of cell IA=1
    double offi;

    // minutes of latitude for non-polar cells on grid A
    double dlat;

//protected:
    blitz::Array<double,1> _dxyp;

protected:
    void init();    // Use after constructor or ncio()

public:
    int size() const { return im * jm; }
    int ndata() const { return size(); }    // Convention makes this more like regular ModelE grids

    double dxyp(int j) const { return _dxyp(j); }

    HntrGrid() {}
    HntrGrid(int _im, int _jm, double _offi, double _dlat);

    HntrGrid(HntrGrid const &other);

    template<class TypeT>
    blitz::Array<TypeT, 2> Array() const
        { return blitz::Array<TypeT,2>(im,jm, blitz::fortranArray); }

    void ncio(ibmisc::NcIO &ncio, std::string const &vname);

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
    /** Initialize overlap data structures, get ready to re-grid.
    TODO: Reference, don't copy, these HntrGrid instances. */
    Hntr(HntrGrid const &_A, HntrGrid const &_B, double _DATMIS);

    Hntr(std::array<HntrGrid const *,2> grids, double _DATMIS = 0)
        : Hntr(*grids[1], *grids[0], _DATMIS) {}


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

    // Default function argument for overlap() template below
    template<typename Typ, bool Val>
    struct IncludeConst
    {
        bool operator()(const Typ& y) const
            { return Val; }
    };

public:

    /** Generates the overlap matrix between two Hntr grids.
    An overlap matrix gives the area of the overlap of gridcells
    (regions on a 2D surface) between two grids.
    Compare to the "original" code in regrid1().
    @param accum Destination for the (sparse) overlap matrix.
    @param R Radius of the sphere, used to calculate area of the grid cells.
    @param includeB Template-type function: bool(int ix) returning True
        if the grid cell from gridB of 1-D index ix is to be included in
        the overlap matrix.
    @see regrid1 */
    template<class AccumT, class IncludeT = IncludeConst<int,true>>
    void overlap(
        AccumT &accum,        // The output (sparse) matrix; 0-based indexing
        double const eq_rad,        // Radius of the Earth
        IncludeT includeB = IncludeT());

};    // class Hntr


template<class AccumT, class IncludeT>
void Hntr::overlap(
    AccumT &accum,        // The output (sparse) matrix; 0-based indexing
    double const eq_rad,        // Radius of the Earth
    IncludeT includeB)
{

    // ------------------
    // Interpolate the A grid onto the B grid
    double const R2 = eq_rad*eq_rad;
    for (int JB=1; JB <= Bgrid.jm; ++JB) {
        int JAMIN = JMIN(JB);
        int JAMAX = JMAX(JB);

        // Buffer for unscaled matrix elements for a single B gridcell
        std::vector<std::pair<int,double>> bvals;

        for (int IB=1; IB <= Bgrid.im; ++IB) {
            int const IJB = IB + Bgrid.im * (JB-1);
            if (!includeB(IJB-1)) continue;

            bvals.clear();
            double WEIGHT = 0;

            int const IAMIN = IMIN(IB);
            int const IAMAX = IMAX(IB);
            for (int JA=JAMIN; JA <= JAMAX; ++JA) {
                double G = SINA(JA) - SINA(JA-1);
                if (JA==JAMIN) G -= GMIN(JB);
                if (JA==JAMAX) G -= GMAX(JB);

                for (int IAREV=IAMIN; IAREV <= IAMAX; ++IAREV) {
                    int const IA  = 1 + ((IAREV-1) % Agrid.im);
                    int const IJA = IA + Agrid.im * (JA-1);

                    double F = 1;
                    if (IAREV==IAMIN) F -= FMIN(IB);
                    if (IAREV==IAMAX) F -= FMAX(IB);

                    double const wt = F*G;
                    WEIGHT += wt;
                    bvals.push_back(std::make_pair(IJA-1, wt));
                }
            }

            // Scale the values we just constructed
            double const byWEIGHT = 1. / WEIGHT;
            for (auto ii=bvals.begin(); ii != bvals.end(); ++ii) {
            	    double const area = R2*Bgrid.dxyp(JB);
                accum.add({IJB-1, ii->first},
                    ii->second * byWEIGHT * area );
            }
        }
    }
}


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

/** Helper for Hntr::overlap() */
class DimClip {
    SparseSetT const *dim;
public:
    DimClip(SparseSetT const *_dim) : dim(_dim) {}

    bool operator()(int ix) const
        { return dim->in_sparse(ix); }
};



}}

#endif    // guard
