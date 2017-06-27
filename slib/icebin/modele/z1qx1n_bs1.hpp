#ifndef ICEBIN_Z1QX1N_BS1_HPP
#define ICEBIN_Z1QX1N_BS1_HPP

#include <cstdio>
#include <boost/algorithm/string.hpp>
#include <icebin/modele/hntr.hpp>


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


extern HntrGrid const g2mx2m;
extern HntrGrid const g10mx10m;
extern HntrGrid const ghxh;
extern HntrGrid const g1x1;
double const dLATM = 180.*60./JM;  //  latitude spacing (minutes)
extern HntrGrid const g1qx1;


#if 0

#define ALL Range::all()

/** Area of memory where a TOPO-generating procedure can place its outputs.
Should be pre-allocated before the generator is called. */
template<class TypeT, int RANK, size_t NVAR>
public ArrayBundle {
    std::string const names;
    std::array<blitz::Array<double, RANK>, NVAR> vars;

    ArrayBundle(
        std::array<std::string, NVAR> const &_names,
        std::array<blitz::Array<double, RANK> *, NVAR> const &_vars,
        std::vector<int> const &shape_v) :
    names(_names)
    {
        auto shape_t(vector_tiny<int,RANK>(shape_v));

        for (size_t i=0; i<NVAR; ++i) {
            if (_vars[i]) {
                // Variable already exists; reference it
                check_dimensions(
                vars[i].reference(*_vars[i]);

                check_dimensions(names[i], vars[i], dims_v);
            } else {
                // Variables does not exist; allocate it!
                vars[i].reference(blitz::Array<TypeT,RANK>(shape_t));
            }
        }
    }

    ArrayBundle<TypeT, RANK, NVAR> c_to_f()
    {
        ArrayBundle<TypeT, RANK, NVAR> ret;
        ret.naems = names;
        for (size_t i=0; i<NVAR; ++i) {
            ret.vars[i].reference(c_to_f(vars[i]));
        }
        return ret;
    }

    ArrayBundle<TypeT, RANK, NVAR> f_to_c()
    {
        ArrayBundle<TypeT, RANK, NVAR> ret;
        ret.naems = names;
        for (size_t i=0; i<NVAR; ++i) {
            ret.vars[i].reference(f_to_c(vars[i]));
        }
        return ret;
    }


};

/** Area of memory where a TOPO-generating procedure can place its outputs.
Should be pre-allocated before the generator is called. */
class TopoOutputs {
    ArrayBundle<double, 2, 14> bundle;
public:
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

    TopoOutputs(ArrayBundle<double, 2, 14> &&_bundle) :
        bundle(std::move(_bundle)),
        FOCEAN(bundle.vars[0]),
        FLAKE(bundle.vars[1]),
        FGRND(bundle.vars[2]),
        FGICE(bundle.vars[3]),
        ZATMO(bundle.vars[4]),
        dZOCEN(bundle.vars[5]),
        dZLAKE(bundle.vars[6]),
        dZGICE(bundle.vars[7]),
        ZSOLDG(bundle.vars[8]),
        ZSGLO(bundle.vars[9]),
        ZLAKE(bundle.vars[10]),
        ZGRND(bundle.vars[11]),
        ZSGHI(bundle.vars[12]),
        FOCENF(bundle.vars[13]) {}
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
    ArrayBundle<double, 2, 8> bundle;
public:
    blitz::Array<double, 2> &FOCEN2, &ZETOP2;
    blitz::Array<double, 2> &FLAKES, &dZGICH, &FGICEH, &ZSOLDH;
    blitz::Array<double, 2> &FCONT1, &FGICE1;

    TopoInputs() :
        FOCEN2(bundle.vars[0]),
        ZETOP2(bundle.vars[1]),
        FLAKES(bundle.vars[2]),
        dZGICH(bundle.vars[3]),
        FGICEH(bundle.vars[4]),
        ZSOLDH(bundle.vars[5]),
        FCONT1(bundle.vars[6]),
        FGICE1(bundle.vars[7])
    {
        // Allocate
        FOCEN2.reference(blitz::Array<double,2>(IM2,JM2,fortranArray));
        ZETOP2.reference(blitz::Array<double,2>(IM2,JM2,fortranArray));
        FLAKES.reference(blitz::Array<double,2>(IMS,JMS,fortranArray));
        dZGICH.reference(blitz::Array<double,2>(IMH,JMH,fortranArray));
        FGICEH.reference(blitz::Array<double,2>(IMH,JMH,fortranArray));
        ZSOLDH.reference(blitz::Array<double,2>(IMH,JMH,fortranArray));
        FCONT1.reference(blitz::Array<double,2>(IM1,JM1,fortranArray));
        FGICE1.reference(blitz::Array<double,2>(IM1,JM1,fortranArray));
    }
};
// =================================================================
/** Read a raw real*4 array and convert to real*8 */
template<int RANK>
void fread_raw4(FILE *finm blitz::Array<double,RANK> &var)
{
    blitz::Array<double,1> var1(
        reshape<double,RANK,1>(var, blitz::shape(-1)));
    blitz::Array<float,1> var4(var.extent(0));
    fread(var4.data(), sizeof(float), nrec, fin);
    for (int i=0; i<nrec; ++i) var1(i) = var4(i);
}

/** Reads and then does a Fortran TRIM() */
std::string fread_string(FILE *fin, size_t len)
{
    char line[len+1];
    fread(line, 1, 80, fin);
    for (size_t i=len-1; i >=0; --i) {
        if (line[i] != ' ') break;
    }
    return std::string(line, i+1);

}

/** Reads and then does a Fortran TRIM() */
std::string fread_string(FILE *fin, char *line, size_t len)
{
    fread(line, 1, 80, fin);
    for (size_t i=len-1; i >=0; --i) {
        if (line[i] != ' ') break;
    }
    line[i+1] = '\0';
}

// https://stackoverflow.com/questions/24611215/one-liner-for-raii-on-non-pointer
struct FILEDeleter
{
    typedef FILE *pointer;
    void operator()(FILE *fp) { fclose(fp); }
};

typedef std::unique_ptr<FILE, FILEDeleter> FilePtr;


void read_raw(TopoInputs &in, std::string const &dir_s)
{
    char titlei[81];    // Extra char for the null character
    boost::filename dir(dir_s);

    printf("BEGIN z1qx1n_bs1 Read Input Files\n");

    // Read in Z2MX2M.NGDC
    {FilePtr fin(fopen((dir / "Z2MX2M.NGDC").string().c_str(), "rb"));
        fread_string(fin, titlei, 80);
        fread_raw4(fin, in.FOCEN2);
        printf("FOCEN2 read from Z2MX2M.NGDC: %s\n", titlei);

        fread_string(fin, titlei, 80);
        fread_raw4(fin, in.ZETOP2);
        printf("ZETOP2 read from Z2MX2M.NGDC: %s\n", titlei);
    }

    // Read in Z10MX10M
    {FilePtr fin(fopen((dir / "Z10MX10M").string().c_str(), "rb"));
        // Read  (2)
        fread_string(fin, titlei, 80);
        fread_raw4(fin, in.FLAKES);
        printf("FLAKES read from Z10MX10M: %s\n", titlei);
    }

    // Read in ZICEHXH
    {FilePtr fin(fopen((dir / "ZICEHXH").string().c_str(), "rb"));
        fread_string(fin, titlei, 80);
        fread_raw4(fin, in.dZGICH);
        printf("dZGICH read from ZICEHXH: %s\n", titlei);

        fread_string(fin, titlei, 80);
        fread_raw4(fin, in.FGICEH);
        printf("FGICEH read from ZICEHXH: %s\n", titlei);

        fread_string(fin, titlei, 80);
        fread_raw4(fin, in.ZSOLDH);
        printf("ZSOLDH read from ZICEHXH: %s\n", titlei);
    }

    // Read in ZNGDC1
    {FilePtr fin(fopen((dir / "ZNGDC1").string().c_str(), "rb"));
        // Read  (4)
        // Read  (4)
        // Read  (4)
        fread_string(fin, titlei, 80);
        fread_raw4(fin, in.FCONT1);
        printf("FCONT1 read from ZNGDC1: %s\n", titlei);

        fread_string(fin, titlei, 80);
        fread_raw4(fin, in.FGICE1);
        printf("FGICE1 read from ZNGDC1: %s\n", titlei);
    }

    printf("END z1qx1n_bs1 Read Input Files\n");
}


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







class ReadZ2Mx2M_NGDC {
    blitz::Array<double, 2> FOCEN2, ZETOP2;
    ReadZ2Mx2M_NGDC() :
        FOCEN2(IM2,JM2),
        ZETOP2(IM2,JM2)
    {
    }
};


blitz::Array<double,2> regrid(
    Hntr &hntr,
    blitz::Array<double,2> const &WTA,
    blitz::Array<double,2> const &A)
{
    blitz::Array<double,2> ret(Bgrid.Array());
    blitz::Array<double,1> ret1(reshape1(ret));
    blitz::Array<double,1> WTA1(reshape1(WTA));
    blitz::Array<double,1> A1(reshape1(A));
    hntr.regrid(WTA1, A1, ret1);
    return ret;
}





void callZ(
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
    blitz::Array<double,2> &ZSGHI)
{
    //
    // Input:  FOCEN2 = ocean fraction at 2 x 2 (minute)
    //         ZSOLD2 = solid topography (above ice) at 2 x 2 (minute)
    //         ZSOLG2 = solid ground topography at 2 x 2 (minute)
    //
    // Output: ZATMO  = atmospheric topography (m)
    //         dZLAKE = mean lake thickness (m)
    //         ZSOLDG = solid ground topography (m)
    //         ZSGLO  = lowest value of ZSOLD2 in model cell (m)
    //         ZLAKE  = surface lake topography (m)
    //         ZGRND  = altitude break between ground and land ice (m)
    //         ZSGHI  = highes value of ZSOLD2 in model cell (m)
    //

    // Local variables
      Integer*4 I,J,IMAX, I2,I21,I2M,J2,J21,J2M, N,NM,NNEW, NLAKE,NGRND
      Real*4    ZN(1296000),AN(1296000)
      Real*8    SAREA,SAZSG, ANSUM,VNSUM, ZLBOT, AZATMO
    //
      If (IM2*JM2/JM > 1296000)  GoTo 800


    for (int J=1; J <= JM; ++J) {
        J21 = (J-1)*JM2/JM + 1;    // 2-minute cells inside (I,J)
        J2M = J*JM2/JM;
        int const IMAX= (J==1 || J==JM ? 1 : IM);
        for (int I=1; I<=IMAX; ++I) {
            int I21 = (I-1)*IM2/IM + 1;
            int I2M = (IMAX == 1 ? IM2 : I*IM2/IM);

            if (FOCEAN(I,J) != 0) {   // (I,J) is an ocean cell
                ZATMO(I,J) = 0;
                dZLAKE(I,J) = 0;
                // ZSOLDG(I,J) = - dZOCEN(I,J)  !  already filled in
                ZLAKE(I,J) = 0;
                ZGRND(I,J) = 0;
                ZSGHI(I,J) = 0;
                ZSGLO(I,J) = 999999;
                for (int J2=j21; J2 <= J2M; ++J2) {
                for (int I2=j21; I2 <= I2M; ++I2) {
                    if (ZSGLO(I,J) > ZSOLD2(I2,J2) && FOCEN2(I2,J2) == 1.) {
                        ZSGLO(I,J) = ZSOLD2(I2,J2);
                    }
                }}

                If (ZSGLO(I,J) == 999999)  GoTo 811
            } else {  // (I,J) is a continental cell
                // Order 2-minute continental cells within (I,J) and sum their area
                struct AreaDepth {
                    double area;
                    double depth;
                    bool operator<(AreaDepth const &other)
                        { return depth < other.depth; }
                    AreaDepth(double _area, double _depth) :
                        area(_area), depth(_depth) {}
                };

                std::vector<AreaDepth> cells2;
                double SAREA = 0;
                double SAZSG = 0;
                int NM = 0;
                for (int J2=j21; J2 <= J2M; ++J2) {
                for (int I2=j21; I2 <= I2M; ++I2) {
                    If (FOCEN2(I2,J2) == 1) continue;
                    double area = g2mx2m.dxyp(J2);
                    cells2.push_back(AreaDepth(area, ZSOLD2(I2,J2)));

                    SAREA += area;
                    SAZSG += area*ZSOLG2(I2,J2)
                }}
                std::sort(cells2.begin(), cells2.end());

                If (SAREA == 0) (*icebin_error)(-1,
                    "Continental cell (%d,%d) has no continental area on 2-minute grid. (%d-%d, %d-%d)",
                    I,J,I21,I2M,J21,J2M);

                // Determine ZSOLDG
                ZSOLDG(I,J) = SAZSG / SAREA
                // Determine ZSGLO
                ZSGLO(I,J) = cells2[0].depth;
                // Determine ZLAKE and dZLAKE
                double ANSUM = 0;  // accumulated area and volume
                double VNSUM = 0;
                int NLAKE;
                for (NLAKE=0; ; ++NLAKE) {
                    if (NLAKE == cells2.size()) {
                        --NLAKE;
                        break;
                    }
                    ANSUM += cells2[NLAKE].area;
                    VNSUM += cells2[NLAKE].area * cells2[nlake].depth;
                    if (ANSUM > SAREA*out.FLAKE(I,J)) break;
                }

                ZLAKE(I,J) = cells2[NLAKE].depth;

                If (out.FLAKE(I,J) > 0) {
                   ZLBOT = (VNSUM - (ANSUM - SAREA*FLAKE(I,J))*cells2[NLAKE].depth) /
                          (SAREA*FLAKE(I,J))
                   dZLAKE(I,J) = Max (cells2[NLAKE].depth-ZLBOT, 1d0)
                } else {
                   dZLAKE(I,J) = 0
                }

                // Determine ZATMO [m]
                AZATMO = ANSUM*ZLAKE(I,J)
                for (int N=NLAKE+1; N<cells2.size(); ++N) {
                    AZATMO += cells2[N].area * cells2[N].depth;
                }
                ZATMO(I,J) = AZATMO / SAREA

                // Determine ZGRND
                ANSUM -= cells2[NLAKE].area;
                int NGRND;
                for (NGRND=NLAKE; ; ++NGRND) {
                    // At end of loop, NGRND points to the last valid item
                    if (NGRND == cells2.size()) {
                        --NGRND;
                        break;
                    }
                    ANSUM += cells2[NGRND].area;
                    If (ANSUM > SAREA*(FLAKE(I,J)+FGRND(I,J))) break;
                }
                ZGRND(I,J) = cells[NGRND].depth;
                // Determine ZSGHI
                ZSGHI(I,J) = cells2[cells.size()-1].depth;
            }
        }
        // Replicate Z data to all longitudes at poles
        If (J==1 || J==JM) {
             ZATMO (Range(2,IM),J) = ZATMO (1,J)
             dZLAKE(Range(2,IM),J) = dZLAKE(1,J)
             ZSOLDG(Range(2,IM),J) = ZSOLDG(1,J)
             ZSGLO (Range(2,IM),J) = ZSGLO (1,J)
             ZLAKE (Range(2,IM),J) = ZLAKE (1,J)
             ZGRND (Range(2,IM),J) = ZGRND (1,J)
             ZSGHI (Range(2,IM),J) = ZSGHI (1,J)
        }
    }
}

void z1qx1n_bs1(TopoInputs &in, TopoOutput &out)
{
    HntrGrid const g2mx2m(HntrGrid(IM2, JM2, 0., 2.));
    HntrGrid const g10mx10m(HntrGrid(IMS, JMS, 0., 10.));
    HntrGrid const ghxh(HntrGrid(IMH, JMH, 0., 30.));
    HntrGrid const g1x1(HntrGrid(IM1, JM1, 0., 60.));
    double const dLATM = 180.*60./JM;  //  latitude spacing (minutes)
    HntrGrid const g1qx1(HntrGrid(IM, JM, 0., dLATM));

    double const TWOPI = 2. * M_PI;
    double const AREAG = 4. * M_PI;

    // Create weight vector of all 1's
    blitz::Array<double, 2> WT2(IM2,JM2);
    WT2 = 1;

    //
    // Add small ice cap and glacier data to FGICEH and dZGICH
    // north of Antarctic area.
    // Continental cells north of 78N are entirely glacial ice.
    FGICE1(ALL, Range(JM1*14/15+1, JM1)) = FCONT1(ALL, Range(JM1*14/15+1, JM1));

    Hntr hntr1h(g1x1m, ghxh);
    auto FCON1H(regrid(hntr1h, WT2, in.FCONT1));
    auto FGIC1H(regrid(hntr1h, WT2, in.FGICE1));

    // RGIC1H = areal ratio of glacial ice to continent
    // For smaller ice caps and glaciers, dZGICH = CONSTK * RGIC1H^.3
    // Constant is chosen so that average value of dZGICH is 264.7 m
    // 264.7  =  sum(DXYP*FGIC1H*dZGICH) / sum(DXYP*FGIC1H)  =
    //        =  CONSTK * sum(DXYP*FGIC1H*RGIC1H^.3) / sum(DXYP*FGIC1H)
    double SUMDFR = 0;
    double SUMDF = 0;
    blitz::Array<double,2> RGIC1H(IMH, JMH);
    for (int JH=1+JMH/6; JH<=JMH; ++JH) {
        for (int IH=1; IH<=IMH; ++IH) {
            if (FGICEH(IH,JH) > 0)
                FGIC1H(IH,JH) = 0;  //  ignore Greenland
            RGIC1H(IH,JH) = FGIC1H(IH,JH) / (FCON1H(IH,JH)+1e-20);
        }
        SUMDFR += ghxh.dxyp(JH) * blitz::sum(FGIC1H(ALL,JH)*RGIC1H(ALL,JH)**.3)
        SUMDF += ghxh.dxyp(JH) * blitz::sum(FGIC1H(ALL,JH))
    }
    double CONSTK = 264.7 * SUMDF / SUMDFR;

    // Replace FGICEH and dZGICH away from Greenland
    for (int JH=JMH/6+1; JH <= JMH; ++JH) {
        for (int IH=1; IH <= IMH; ++IH) {
            If (FGICEH(IH,JH) == 0) {
                FGICEH(IH,JH) = FGIC1H(IH,JH)
                dZGICH(IH,JH) = CONSTK * RGIC1H(IH,JH)**.3
            }
        }
    }


    // ETOPO2 treats Antarctic ice shelves as ocean.
    // When this happens ETOPO2 data are replaced with interpolated data
    // from FGICEH and dZGICH.  Resulting files are:
    // FOCEN2 = Ocean fraction (0 or 1) correct for Antarctic ice shelves
    // FCONT2 = Continental fraction (0 or 1)
    // FGICE2 = Glacial ice fraction (0 or 1)
    // dZGIC2 = Thickness of glacial ice (m)
    // ZSOLD2 = Solid topography (m)        (above ice)
    // ZSOLG2 = Solid ground topography (m) (beneath ice)
    //
    Hntr hntrhm2(ghxh, g2mxg2m);
    auto FGICE2(regrid(hntrhm2, WT2, FGICEH));    // WT2 is too big...
    auto dZGIC2(regrid(hntrhm2, FGICEH, dZGICH));
    auto ZSOLD2(regrid(hntrhm2, FGICEH, ZSOLDH));

    // North of Antarctic area: 60S to 90N
    blitz::Array<double,2> FCONT2(IM2, JM2);
    blitz::Array<double,2> ZSOLG2(IM2, JM2);
    for (int J2=JM2/6+1; J2 <= JM2; ++J2) {
        FCONT2(ALL, J2) = 1. - in.FOCEN2(ALL, J2);
        FGICE2(ALL, J2) = FGICE2(ALL, J2) * FCONT2(ALL, J2);
        dZGIC2(ALL, J2) = dZGIC2(ALL, J2) * FCONT2(ALL, J2);
        ZSOLD2(ALL, J2) = ZETOP2(ALL, J2);
        ZSOLG2(ALL, J2) = ZETOP2(ALL, J2) - dZGIC2(ALL, J2);
    }

    // Antarctic area: 90S to 60S
    for (int J2=1; J2<=JM2/6; ++J2) {
    for (int I2=1; I2<=JM2; ++I2) {
        if (in.FOCEN2(I2,J2) == 0) {
            // in.FOCEN2(I2,J2) = 0
            FCONT2(I2,J2) = 1;
            FGICE2(I2,J2) = 1;
            // dZGIC2(I2,J2) = dZGIC2(I2,J2)  //  ZETOP2 has 2m and other low
            // ZSOLD2(I2,J2) = ZSOLD2(I2,J2)       values over ice shelves
            if (ZETOP2(I2,J2) >= 100)  ZSOLD2(I2,J2) = ZETOP2(I2,J2);
            ZSOLG2(I2,J2) = ZSOLD2(I2,J2) - dZGIC2(I2,J2);
        } else if (FGICE2(I2,J2) <= .5) {  //  and in.FOCEN2(I2,J2) == 1
            // in.FOCEN2(I2,J2) = 1
            FCONT2(I2,J2) = 0;
            FGICE2(I2,J2) = 0;
            dZGIC2(I2,J2) = 0;
            ZSOLD2(I2,J2) = ZETOP2(I2,J2);
            ZSOLG2(I2,J2) = ZETOP2(I2,J2);
        } else if (FGICE2(I2,J2) > .5) {  //  and in.FOCEN2(I2,J2) == 1
            in.FOCEN2(I2,J2) = 0;
            FCONT2(I2,J2) = 1;
            FGICE2(I2,J2) = 1;
            dZGIC2(I2,J2) = ZSOLD2(I2,J2) - ZETOP2(I2,J2);
            // ZSOLD2(I2,J2) = ZSOLD2(I2,J2)
            ZSOLG2(I2,J2) = ZETOP2(I2,J2);
        }
    }}


    //
    // FOCEAN: Ocean Surface Fraction (0:1)
    //
    // Fractional ocean cover FOCENF is interpolated from FOAAH2
    Hntr hntr2mq1(g2mx2m, g1qx1);
    auto FOCENF(regridp(hntr2mq1, WT2, FOCEN2));    // Fractional ocean cover

    // FOCEAN (0 or 1) is rounded from FOCEAN
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        FOCEAN(i,j) = std::round(FOCENF(i,j));
    }}

    // Set following grid cells to be continent
    FOCEAN( 84, 18) = 0;
    FOCEAN( 85, 18) = 0;
    FOCEAN(236, 82) = 0;
    FOCEAN(242, 82) = 0;
    FOCEAN(245, 82) = 0;
    FOCEAN(224,101) = 0;
    FOCEAN( 53,119) = 0;
    FOCEAN(171,125) = 0;  //  Cyprus
    FOCEAN(164,126) = 0;  //  Crete
    FOCEAN(158,129) = 0;
    FOCEAN(158,130) = 0;
    FOCEAN(242,131) = 0;
    FOCEAN(263,136) = 0;
    FOCEAN(258,137) = 0;
    FOCEAN(258,138) = 0;
    FOCEAN( 46,139) = 0;
    FOCEAN(258,139) = 0;
    FOCEAN(275,152) = 0;
    FOCEAN(  8,156) = 0;
    FOCEAN( 10,156) = 0;
    FOCEAN( 12,157) = 0;
    FOCEAN(172,157) = 0;
    FOCEAN(202,157) = 0;
    FOCEAN( 69,159) = 0;
    FOCEAN(204,159) = 0;
    FOCEAN( 62,167) = 0;
    FOCEAN( 73,171) = 0;
    FOCEAN( 75,171) = 0;
    FOCEAN( 78,171) = 0;

    // Set following grid cells to be ocean
    FOCEAN(179,105) = 1;
    FOCEAN( 54,119) = 1;
    FOCEAN(241,131) = 1;
    FOCEAN(258,143) = 1;
    FOCEAN(165,150) = 1;
    FOCEAN(274,152) = 1;
    FOCEAN( 15,154) = 1;
    FOCEAN( 92,155) = 1;
    FOCEAN( 13,157) = 1;
    FOCEAN(173,157) = 1;
    FOCEAN(176,157) = 1;
    FOCEAN(203,157) = 1;
    FOCEAN( 55,159) = 1;
    FOCEAN(103,159) = 1;
    FOCEAN(203,159) = 1;
    FOCEAN( 67,160) = 1;
    FOCEAN( 68,160) = 1;
    FOCEAN( 79,160) = 1;
    FOCEAN(199,160) = 1;
    FOCEAN(126,161) = 1;
    FOCEAN( 68,162) = 1;
    FOCEAN( 75,165) = 1;
    FOCEAN(225,169) = 1;

    // Average non-fractional and fractional ocean covers over latitude
    printf(
        " Comparison between Fractional and Non-fractional Ocean Cover\n\n"
        "         # of      # of     differ\n"
        "         fract    NOfrac      in #\n"
        "   J     cells     cells     cells\n"
        "   =     =====     =====     =====\n");
    blitz::Array<double,2> FOFLAT(JM);
    blitz::Array<double,2> FONLAT(JM);
    for (int J=JM; J >= 1; --J) {
        FOFLAT(J) = blitz::sum(FOCENF(ALL,J));
        FONLAT(J) = blitz::sum(FOCEAN(ALL,J));
    }
    double const byAREAG = 1. / (4.*M_PI);
    printf("%4d%102f%10.2f%10.2f\n", J,FOFLAT(J),FONLAT(J),FOFLAT(J)-FONLAT(J));
    double FOFSH = IM*JM *
        blitz::sum(FOFLAT(Range(1,JM/2)) * g1qx1.dxyp(Range(1,JM/2))) * byAREAG;
    double FONSH = IM*JM *
        blitz::sum(FONLAT(Range(1,JM/2)) * g1qx1.dxyp(Range(1,JM/2))) * byAREAG;
    double FOFNH = IM*JM *
        blitz::sum(FOFLAT(Range(1+JM/2,JM)) * g1qx1.dxyp(Range(1+JM/2,JM))) * byAREAG;
    double FONNH = IM*JM *
        blitz::sum(FONLAT(Range(1+JM/2,JM)) * g1qx1.dxyp(Range(1+JM/2,JM))) * byAREAG;
    printf("NH: %f %f %f\n", FOFNH, FONNH, FOFNH-FONNH);
    printf("SH: %f %f %f\n", FOFSH, FONSH, FOFSH-FONSH);


    //
    // FLAKE: Lake Surface Fraction (0:1)
    //
    // FLAKE is interpolated from FLAKES
    Hntr hntr10m1q(g10mx10m, g1qx1);
    regrid1(hntr10m1q, WT2, in.FLAKES, out.FLAKE);

    // Antarctica and Arctic area have no lakes
    out.FLAKE(ALL, 1:JM/6) = 0;             //  90:60 S
    out.FLAKE(ALL, JM*14/15+1:JM) = 0;      //  78:90 N
    out.FLAKE(:IM/2,JM*41/45+1:JM) = 0;  //  74:90 N, 0:180 W

    for (int J=JM*5/6 J <= JM*11/12; ++J) {    //  southern
    for (int I=IM/3+1; I <= (int)(.5 + .75*IM*(J-JM*.3)/JM); ++I) {  //  Greenland
       out.FLAKE(I,J) = 0;
    }}

    // Apportion out.FLAKE to the nonocean fraction and round to 1/256
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        FLAKE(i,j) = FLAKE(i,j)*(1-FOCEAN(i,j)) / (1-FOCENF(i,j)+1e-20);
        FLAKE(i,j) = std::round(FLAKE(i,j)*256) / 256.;
    }}


    //
    // FGICE: Glacial Ice Surface Fraction (0:1)
    //
    // FGICE is interpolated from FGICE2
    regrid(hntr2mq1, FCONT2, FGICE2, out.FGICE);

    // Antarctica is entirely glacial ice, no lakes nor ground
    out.FGICE(ALL,Range(1,JM/6)) = 1. - out.FOCEAN(ALL,Range(1,JM/6));

    // Continental cells north of 78N are entirely glacial ice
    out.FGICE(ALL,Range(JM*14/15+1,JM)) = 1 - out.FOCEAN(ALL,Range(JM*14/15+1,JM));

    // There is no glacial ice over oceans
    out.FGICE(ALL,Range(JM/6+1,JM)) =
        out.FGICE(ALL,Range(JM/6+1,JM)) * (1-out.FOCEAN(ALL,Range(JM/6+1,JM)));

    // Round out.FGICE to nearest 1/256
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        out.FGICE(i,j) = std::round(out.FGICE(i,j)*256) / 256.;
    }}

    // Check that out.FGICE is between 0 and 1
    // If out.FGICE+FLAKE exceeds 1, reduce FLAKE
    for (int J=JM/6+1; J <= JM; ++J) {
    for (int I=1; I<=IM; ++I) {
        if (out.FGICE(I,J) < 0) {
            printf("210: out.FGICE(%d,%d) < 0: %g\n" ,I,J,out.FGICE(I,J));
            out.FGICE(I,J) = 0;
        }
        if (out.FGICE(I,J) > 1) {
            printf("210: out.FGICE(%d,%d) > 1: %g\n" ,I,J,out.FGICE(I,J));
            out.FGICE(I,J) = 1;
        }
        if (out.FLAKE(I,J)+out.FGICE(I,J)+out.FOCEAN(I,J) > 1) {
            printf("210: FGICE+FLAKE+FOCEAN (%d,%d) > 1: %g + %g + %g\n",
                I,J,out.FGICE(I,J), out.FLAKE(I,J),out.FOCEAN(I,J));
            out.FLAKE(I,J) = 1. - out.FGICE(I,J) - out.FOCEAN(I,J);
        }
    }}

    // Replace land cells without vegetation with glacial ice in Z1QX1N
    // out.FGICE(35,52) = 1 - FLAKE(35,52)     4x3 Model

    //
    // FGRND: Surface Fraction of Ground (0:1)
    //
    out.FGRND = 1 - out.FOCEAN - out.FLAKE - out.FGICE;

    // Check that FGRND is between 0 and 1
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        if (out.FGRND(i,j) < 0 || out.FGRND(i,j) > 1) {
            printf("Error: FGRND(%d,%d) = %g %g %g %g\n", i,j,
                out.FGRND(i,j), out.FOCEAN(i,j), out. FLAKE(i,j), out.FGICE(i,j));
        }
    }}

    //
    // dZOCEN: Ocean Thickness (m)
    //
    auto dZOCEN(regrid(hntr2mq1, FOCEN2, ZSOLG2));
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        dZOCEN(i,j) = -dZOCEN(i,j) * out.FOCEAN(i,j);

        // Check that dZOCEN is positive
        If (out.FOCEAN(i,j) == 1 && dZOCEN(i,j) <= 0) {
            printf("Error: dZOCEN(%d,%d) <= 0 %g",i,j,dZOCEN(i,j));
        }
    }}

    //
    // dZGICE: Glacial Ice Thickness (m)
    //
    auto dZGICE(regrid(hntr2mq1, FGICE2, dZGICE2));
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        If (out.FGICE(I,J) > 0) {
            dZGICE(I,J) = std::max(dZGICE(I,J), 1.);
        } else {
            dZGICE(I,J) = 0;
        }
    }}

    //
    // ZATMO  = Atmospheric topography (m)
    // dZLAKE = Mean lake thickness (m)
    // ZSOLDG = Solid ground topography (m)
    // ZSGLO  = Lowest value of ZSOLD2 in model cell (m)
    // ZLAKE  = Surface lake topography (m)
    // ZGRND  = Altitude break between ground and land ice (m)
    // ZSGHI  = Highest value of ZSOLD2 in model cell (m)
    //
    out.ZSOLDG = - dZOCEN;  //  solid ground topography of ocean
    callZ(
        FOCEN2,ZSOLD2,ZSOLG2,
        out.FOCEAN, out.FLAKE, out.FGRND, out.ZATMO,
        dZLAKE,ZSOLDG,ZSGLO,ZLAKE,ZGRND,ZSGHI);

    // Reset ZATMO, dZOCEN and ZLAKE by hand if FLAKE(I,J) == 1
    //

    static std::vector<std::tuple<double,    // elev
        std::vector<std::array<int,2>>>> const resets
    {
        // Caspian Sea
        {-30., {    // Elevation [m]
            {186,128}, {187,128}, {185,129}, {186,129},
            {185,130}, {186,130}, {186,131}, {185,132},
            {184,133}, {185,133}, {184,134}, {183,135}, {184,135}
        }},

        // Aral Sea
        {53., {    // Elevation [m]
            {192,135}
        }},
        // Lake Superior
        {75., {    // Elevation [m]
            {75,138}
        }}
    };

    for (auto &reset : rests) {
        double elev = reset.get<0>();
        auto &ijs(reset.get<1>());

        for (std::array<int,2> const &ij : ijs) {
            int const i(ij[0]);
            int const j(ij[1]);
            dZLAKE(i,j) += (elev - ZLAKE(i,j));
            ZATMO(i,j) = elev;
            ZLAKE(i,j) = elev;
        }
    }

    for (int J=1; J <= JM; ++J) {
    for (int I=1; I<=IM; ++I) {
        If (FLAKE(I,J) == 1) {
            printf("FLAKE(%d,%d) == 1: %g %g %g\n",
                I, J, ZATMO(I,J),dZLAKE(I,J),ZLAKE(I,J));
        }
    }}
}
#endif

}}
#endif     // guard
