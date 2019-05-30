#ifndef ICEBIN_HNTR_HPP
#define ICEBIN_HNTR_HPP

#include <ibmisc/blitz.hpp>
#include <ibmisc/indexing.hpp>
#include <icebin/eigen_types.hpp>
#include <icebin/GridSpec.hpp>

namespace icebin {
namespace modele {

extern blitz::Array<double,1> make_dxyp(
    HntrSpec const &spec,
    blitz::GeneralArrayStorage<1> const &storage = blitz::GeneralArrayStorage<1>());


class HntrGrid {
public:
    HntrSpec spec;
//    int im;    // Number of cells in east-west direction
//    int jm;    // Number of cells in north-south direction

//    // number (fraction) of cells in east-west direction from
//    // International Date Line (180) to western edge of cell IA=1
//    double offi;

//    // minutes of latitude for non-polar cells on grid A
//    double dlat;

protected:
    void init();    // Use after constructor or ncio()

public:
    /** Area of grid cell at lattitude index = j.
    @param j One-based indexing. */
    blitz::Array<double,1> dxyp;

    /** Standardized description of indexing for ModelE grids */
    ibmisc::Indexing indexing;

#if 0
    int size() const { return spec.size(); }
    int ndata() const { return size(); }    // Convention makes this more like regular ModelE grids
#endif

    HntrGrid() {}
    HntrGrid(int _im, int _jm, double _offi, double _dlat);

    HntrGrid(HntrGrid const &other);
    explicit HntrGrid(HntrSpec const &spec);

    void ncio(ibmisc::NcIO &ncio, std::string const &vname);

};

template<class TypeT>
blitz::Array<TypeT,2> hntr_array(HntrSpec const &spec)
    { return blitz::Array<TypeT,2>(spec.im,spec.jm, blitz::fortranArray); }


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
//    Hntr(HntrGrid const &_A, HntrGrid const &_B, double _DATMIS);

//    Hntr(std::array<HntrGrid const *,2> grids, double _DATMIS = 0)
//        : Hntr(*grids[1], *grids[0], _DATMIS) {}

    Hntr(double yp17, HntrSpec const &_B, HntrSpec const &_A, double _DATMIS=0.0);




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
    template<class WeightT, class SrcT, class DestT, int RANK>
    void regrid(
        blitz::Array<WeightT,RANK> const &_WTA,
        blitz::Array<SrcT,RANK> const &_A,
        blitz::Array<DestT,RANK> const &B,
        bool mean_polar=false,
        double wtm=1.0, double wtb=0.0) const;


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


    /** Generalized regridding "engine."  used to implement overlap()
        and scaled_regrid_matrix() */
    template<class MatAccumT, class IncludeT>
    void matrix(
        MatAccumT &&mataccum,        // The output (sparse) matrix; 0-based indexing
        IncludeT includeB) const;

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
        AccumT &&accum,        // The output (sparse) matrix; 0-based indexing
        double const eq_rad,        // Radius of the Earth
        IncludeT includeB = IncludeT());

    /** Produces a scaled regrid matrix, without the extra baggage.
    Equivalent to running overlap() and then scaling. */
    template<class AccumT, class IncludeT = IncludeConst<int,true>>
    void scaled_regrid_matrix(
        AccumT &&accum,        // The output (sparse) matrix; 0-based indexing
        IncludeT includeB = IncludeT());
};    // class Hntr



// --------------------------------------------------------------

/** Helper for Hntr::overlap() */
class DimClip {
    SparseSetT const *dim;
public:
    DimClip(SparseSetT const *_dim) : dim(_dim) {}

    bool operator()(int ix) const
        { return dim->in_sparse(ix); }
};


// ----------------------------------------------------------------

// ==================================================================
template<class MatAccumT, class IncludeT>
void Hntr::matrix(
    MatAccumT &&mataccum,        // The output (sparse) matrix; 0-based indexing
    IncludeT includeB) const
{
    // ------------------
    // Interpolate the A grid onto the B grid
    for (int JB=1; JB <= Bgrid.spec.jm; ++JB) {
        int JAMIN = JMIN(JB);
        int JAMAX = JMAX(JB);


        for (int IB=1; IB <= Bgrid.spec.im; ++IB) {
            int const IJB = IB + Bgrid.spec.im * (JB-1);
            if (!includeB(IJB-1)) continue;

            mataccum.clear();

            int const IAMIN = IMIN(IB);
            int const IAMAX = IMAX(IB);
            for (int JA=JAMIN; JA <= JAMAX; ++JA) {
                double G = SINA(JA) - SINA(JA-1);
                if (JA==JAMIN) G -= GMIN(JB);
                if (JA==JAMAX) G -= GMAX(JB);

                for (int IAREV=IAMIN; IAREV <= IAMAX; ++IAREV) {
                    int const IA  = 1 + ((IAREV-1) % Agrid.spec.im);
                    int const IJA = IA + Agrid.spec.im * (JA-1);

                    double F = 1;
                    if (IAREV==IAMIN) F -= FMIN(IB);
                    if (IAREV==IAMAX) F -= FMAX(IB);

                    mataccum.addA(IJA, F*G);
                }
            }

            mataccum.finishB(IJB, JB);
        }
    }
}


// ----------------------------------------------------------
template<class AccumT>
class OverlapMatAccum {
    AccumT &accum;
    HntrGrid const &Bgrid;
    double const R2;

    // Buffer for unscaled matrix elements for a single B gridcell
    std::vector<std::pair<int,double>> bvals;
    double WEIGHT;

public:
    OverlapMatAccum(AccumT &&_accum, HntrGrid const &_Bgrid, double _R2) : accum(_accum), Bgrid(_Bgrid), R2(_R2) {}

    void clear()
    {
        bvals.clear();
        WEIGHT = 0.;
    }

    void addA(int const IJA, double const FG)
    {
        WEIGHT += FG;
        bvals.push_back(std::make_pair(IJA-1, FG));
    }

    void finishB(int const IJB, int const JB)
    {
        // Scale the values we just constructed
        double const byWEIGHT = 1. / WEIGHT;
        for (auto ii=bvals.begin(); ii != bvals.end(); ++ii) {
            double const area = R2*Bgrid.dxyp(JB);
            double const val = ii->second * byWEIGHT * area;
            accum.add({IJB-1, ii->first}, val);
        }
    }

};

template<class AccumT, class IncludeT>
void Hntr::overlap(
    AccumT &&accum,        // The output (sparse) matrix; 0-based indexing
    double const eq_rad,        // Radius of the Earth
    IncludeT includeB)
{
    matrix(OverlapMatAccum<AccumT>(std::move(accum), Bgrid, eq_rad*eq_rad), includeB);
}
// ----------------------------------------------------------
template<class AccumT>
class ScaledRegridMatAccum {
    AccumT accum;
    HntrGrid const &Agrid;

    // Buffer for unscaled matrix elements for a single B gridcell
    std::vector<std::pair<int,double>> bvals;
    double WEIGHT;

public:
    ScaledRegridMatAccum(AccumT &&_accum, HntrGrid const &_Agrid) : accum(std::move(_accum)), Agrid(_Agrid) {}

    void clear()
    {
        bvals.clear();
        WEIGHT = 0.;
    }

    void addA(int const IJA, double const FG)
    {
        WEIGHT += FG;
        bvals.push_back(std::make_pair(IJA-1, FG));
    }

    void finishB(int const IJB, int const JB)
    {
        // Scale the values we just constructed
        double const byWEIGHT = 1. / WEIGHT;
        for (auto ii=bvals.begin(); ii != bvals.end(); ++ii) {
            double const val = ii->second * byWEIGHT;
            accum.add({IJB-1, ii->first}, val);
        }
    }

};

template<class AccumT, class IncludeT>
void Hntr::scaled_regrid_matrix(
    AccumT &&accum,        // The output (sparse) matrix; 0-based indexing
    IncludeT includeB)
{
    matrix(
        ScaledRegridMatAccum<AccumT>(std::move(accum), Agrid),
        std::move(includeB));
}
// ---------------------------------------------------
// ----------------------------------------------------------
template<class WeightT, class SrcT, class DestT>
class RegridAccum {
    blitz::Array<WeightT,1> &WTA;
    blitz::Array<SrcT,1> &A;
    blitz::Array<DestT,1> &B;
    double const DATMIS;
    double const wtm;
    double const wtb;

    double WEIGHT;
    double VALUE;

public:
    RegridAccum(
        blitz::Array<WeightT,1> &_WTA,
        blitz::Array<SrcT,1> &_A,
        blitz::Array<DestT,1> &_B,
        double _DATMIS,
        double _wtm, double _wtb)    // Use wtm * WTA + wtb for weight
    : WTA(_WTA), A(_A), B(_B), DATMIS(_DATMIS), wtm(_wtm), wtb(_wtb) {}

    void clear()
    {
        VALUE = 0.;
        WEIGHT = 0.;
    }

    void addA(int const IJA, double const FG)
    {
        double const wta = wtm * WTA(IJA) + wtb;
        double const wt = FG * wta;
        WEIGHT += wt;
        VALUE += wt * A(IJA);
    }

    void finishB(int const IJB, int const JB)
    {
        B(IJB) = (WEIGHT == 0 ? DATMIS : VALUE / WEIGHT);
    }

};


template<class WeightT, class SrcT, class DestT, int RANK>
void Hntr::regrid(
    blitz::Array<WeightT,RANK> const &_WTA,
    blitz::Array<SrcT,RANK> const &_A,
    blitz::Array<DestT,RANK> const &_B,
    bool mean_polar,
    double wtm, double wtb) const
{
    // Reshape to 1-D
    auto WTA(ibmisc::reshape1(_WTA, 1));
    auto A(ibmisc::reshape1(_A, 1));
    auto B(ibmisc::reshape1(_B, 1));

    // Check array dimensions
    if ((WTA.extent(0) != Agrid.spec.size()) ||
        (A.extent(0) != Agrid.spec.size()) ||
        (B.extent(0) != Bgrid.spec.size()))
    {
        (*icebin_error)(-1, "Error in dimensions: (%d, %d, %d) vs. (%d, %d)\n",
            WTA.extent(0), A.extent(0), B.extent(0),
            Agrid.spec.size(), Bgrid.spec.size());
    }


    matrix(
        RegridAccum<WeightT,SrcT,DestT>(WTA, A, B, DATMIS, wtm, wtb),
        IncludeConst<int,true>());

    if (mean_polar) {
        // Replace individual values near the poles by longitudinal mean
        for (int JB=1; JB <= Bgrid.spec.jm; JB += Bgrid.spec.jm-1) {
            double BMEAN  = DATMIS;
            double WEIGHT = 0;
            double VALUE  = 0;
            for (int IB=1; ; ++IB) {
                if (IB > Bgrid.spec.im) {
                    if (WEIGHT != 0) BMEAN = VALUE / WEIGHT;
                    break;
                }
                int IJB = IB + Bgrid.spec.im * (JB-1);
                if (B(IJB) == DATMIS) break;
                WEIGHT += 1;
                VALUE  += B(IJB);
            }
            for (int IB=1; IB <= Bgrid.spec.im; ++IB) {
                int IJB = IB + Bgrid.spec.im * (JB-1);
                B(IJB) = BMEAN;
            }
        }
    }
}

/** Creates a HntrSpec for the Atmosphere grid by halving a HntrSpec
for the Ocean grid.  Relies on this 2-to-1 relationship of ocean to
atmosphere in ModelE. */
extern HntrSpec make_hntrA(HntrSpec const &hntrO);

}}    // namespace
#endif    // guard
