#include <cmath>
#include <ibmisc/fortranio.hpp>
#include <ibmisc/ncbulk.hpp>
#include <icebin/error.hpp>
#include <icebin/modele/make_topoo.hpp>
#include <icebin/modele/grids.hpp>

using namespace blitz;
using namespace ibmisc;
using namespace spsparse;

namespace icebin {
namespace modele {

static double const NaN = std::numeric_limits<double>::quiet_NaN();

/** Within each 1-minute gridcell, it has a hi-res topography and a
lake fraction.  It determines the elevation of the lake, being at the
bottom of the gridcell. */
static void callZ(
    // (IM1m, JM1m)
    blitz::Array<int16_t,2> const &FGICE1m,
    blitz::Array<int16_t,2> const &FOCEAN1m,
    blitz::Array<int16_t,2> const &ZICETOP1m,
    blitz::Array<int16_t,2> const &ZSOLG1m,

    // (IM, IM)
    blitz::Array<double,2> &FOCEAN,
    blitz::Array<double,2> &FLAKE,
    blitz::Array<double,2> &FGRND,

    // (IM, IM)
    blitz::Array<double,2> &ZATMO,
    blitz::Array<double,2> &ZATMOF,
    blitz::Array<double,2> &dZLAKE,
    blitz::Array<double,2> &ZSOLDG,
    blitz::Array<double,2> &ZICETOP,
    blitz::Array<double,2> &ZSGLO,
    blitz::Array<double,2> &ZLAKE,
    blitz::Array<double,2> &ZGRND,
    blitz::Array<double,2> &ZSGHI)
{

    //
    // Input:  FOCEAN1m = ocean fraction at 2 x 2 (minute)
    //         ZICETOP1m = solid topography (above ice) at 2 x 2 (minute)
    //         ZSOLG1m = solid ground topography at 2 x 2 (minute)
    //
    // Output: ZATMO  = atmospheric topography (m)
    //         dZLAKE = mean lake thickness (m)
    //         ZSOLDG = solid ground topography (m)
    //         ZICETOP = atmospheric topography of JUST ice-covered region (m)
    //         ZSGLO  = lowest value of ZICETOP1m in model cell (m)
    //         ZLAKE  = surface lake topography (m)
    //         ZGRND  = altitude break between ground and land ice (m)
    //         ZSGHI  = highes value of ZICETOP1m in model cell (m)
    //

    HntrGrid grid_g1mx1m(g1mx1m);

    for (int J=1; J <= JM; ++J) {
        int J11 = (J-1)*JM1m/JM + 1;    // 1-minute cells inside (I,J)
        int J1M = J*JM1m/JM;
        int const IMAX= (J==1 || J==JM ? 1 : IM);
        for (int I=1; I<=IMAX; ++I) {
            int I11 = (I-1)*IM1m/IM + 1;
            int I1M = (IMAX == 1 ? IM1m : I*IM1m/IM);

            if (FOCEAN(I,J) != 0) {   // (I,J) is an ocean cell
                ZATMO(I,J) = 0;
                dZLAKE(I,J) = 0;
                // ZSOLDG(I,J) = - dZOCEN(I,J)  //  already filled in
                ZLAKE(I,J) = 0;
                ZGRND(I,J) = 0;
                ZSGHI(I,J) = 0;
                ZSGLO(I,J) = 999999;
                ZICETOP(I,J) = 0;
                for (int J2=J11; J2 <= J1M; ++J2) {
                for (int I2=I11; I2 <= I1M; ++I2) {
                    if (ZSGLO(I,J) > ZICETOP1m(I2,J2) && FOCEAN1m(I2,J2) == 1) {
                        ZSGLO(I,J) = ZICETOP1m(I2,J2);
                    }
                }}

                if (ZSGLO(I,J) == 999999) (*icebin_error)(-1,
                    "Ocean cell FOCEAN(%d,%d)=%g has no ocean area on 2-minute grid; "
                    "I11,I1M,J11,J1M = %d,%d,%d,%d",
                    I,J,FOCEAN(I,J),
                    I11,I1M,J11,J1M);

                // --------- Determine ZICETOP and ZATMOF, even for ocean-only cell
                double SAREA = 0;    // Entire gridcell...
                double SAZSG = 0;
                double SAREA_li = 0;    // Landice portions only
                double SAZSG_li = 0;
                for (int J1m=J11; J1m <= J1M; ++J1m) {
                for (int I1m=I11; I1m <= I1M; ++I1m) {
                    double area = grid_g1mx1m.dxyp(J1m);
                    SAREA += area;      // Surface area of LAND (non-ocean)

                    if (FOCEAN1m(I1m,J1m) == 1) continue;
                    SAZSG += area*ZSOLG1m(I1m,J1m);

                    if (FGICE1m(I1m,J1m) != 0) {
                        SAREA_li += area;
                        SAZSG_li += area*ZICETOP1m(I1m,J1m);
                    }
                }}

                ZICETOP(I,J) = (SAREA_li > 0 ? SAZSG_li / SAREA_li : 0);
                ZATMOF(I,J) = SAZSG / SAREA;
            } else {  // (I,J) is acontinent cell
                // Order 1-minute continental cells within (I,J) and sum their area
                struct AreaDepth {
                    double area;
                    double depth;
                    bool operator<(AreaDepth const &other)
                        { return depth < other.depth; }
                    AreaDepth(double _area, double _depth) :
                        area(_area), depth(_depth) {}
                };

                std::vector<AreaDepth> cells2;
                double SAREA = 0;    // Entire gridcell...
                double SAZSG = 0;
                double SAREA_li = 0;    // Landice portions only
                double SAZSG_li = 0;
                int NM = 0;
                for (int J1m=J11; J1m <= J1M; ++J1m) {
                for (int I1m=I11; I1m <= I1M; ++I1m) {
                    if (FOCEAN1m(I1m,J1m) == 1) continue;
                    double area = grid_g1mx1m.dxyp(J1m);
                    cells2.push_back(AreaDepth(area, ZICETOP1m(I1m,J1m)));

                    SAREA += area;      // Surface area of LAND (non-ocean)
                    SAZSG += area*ZSOLG1m(I1m,J1m);

                    if (FGICE1m(I1m,J1m) != 0) {
                        SAREA_li += area;
                        SAZSG_li += area*ZICETOP1m(I1m,J1m);
                    }
                }}
                std::sort(cells2.begin(), cells2.end());

                if (SAREA == 0) (*icebin_error)(-1,
                    "Continental cell (%d,%d) has no continental area on 2-minute grid. (%d-%d, %d-%d)",
                    I,J,I11,I1M,J11,J1M);

                // Determine ZSOLDG
                ZSOLDG(I,J) = SAZSG / SAREA;
                if (SAREA_li > 0) {
                    ZICETOP(I,J) = SAZSG_li / SAREA_li;
                } else {
                    ZICETOP(I,J) = 0;
                }

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
                    VNSUM += cells2[NLAKE].area * cells2[NLAKE].depth;
                    if (ANSUM > SAREA*FLAKE(I,J)) break;
                }

                ZLAKE(I,J) = cells2[NLAKE].depth;

                if (FLAKE(I,J) > 0) {
                   double ZLBOT = (VNSUM - (ANSUM - SAREA*FLAKE(I,J))*cells2[NLAKE].depth) /
                          (SAREA*FLAKE(I,J));
                   dZLAKE(I,J) = std::max(cells2[NLAKE].depth-ZLBOT, 1.0);
                } else {
                   dZLAKE(I,J) = 0;
                }

                // Determine ZATMO [m]
                double AZATMO = ANSUM*ZLAKE(I,J);
                for (int N=NLAKE+1; N<cells2.size(); ++N) {
                    AZATMO += cells2[N].area * cells2[N].depth;
                }
                ZATMO(I,J) = AZATMO / SAREA;

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
                    if (ANSUM > SAREA*(FLAKE(I,J)+FGRND(I,J))) break;
                }
                ZGRND(I,J) = cells2[NGRND].depth;
                // Determine ZSGHI
                ZSGHI(I,J) = cells2[cells2.size()-1].depth;
            }
        }

        // Replicate Z data to all longitudes at poles
        if (J==1 || J==JM) {
             ZATMO (Range(2,IM),J) = ZATMO (1,J);
             dZLAKE(Range(2,IM),J) = dZLAKE(1,J);
             ZSOLDG(Range(2,IM),J) = ZSOLDG(1,J);
             ZICETOP(Range(2,IM),J) = ZICETOP(1,J);
             ZSGLO (Range(2,IM),J) = ZSGLO (1,J);
             ZLAKE (Range(2,IM),J) = ZLAKE (1,J);
             ZGRND (Range(2,IM),J) = ZGRND (1,J);
             ZSGHI (Range(2,IM),J) = ZSGHI (1,J);
        }
    }
}

struct ElevPoints {
    double elev;
    std::vector<std::array<int,2>> points;

    ElevPoints(double _elev, std::vector<std::array<int,2>> &&_points)
        : elev(_elev), points(std::move(_points)) {}
};

static const std::vector<ElevPoints> resets
{
    // Caspian Sea
    ElevPoints(-30., {    // Elevation [m]
        {186,128}, {187,128}, {185,129}, {186,129},
        {185,130}, {186,130}, {186,131}, {185,132},
        {184,133}, {185,133}, {184,134}, {183,135}, {184,135}
    }),

    // Aral Sea
    ElevPoints(53., {    // Elevation [m]
        {192,135}
    }),
    // Lake Superior
    ElevPoints(75., {    // Elevation [m]
        {75,138}
    })
};



static ibmisc::ArrayBundle<double,2> _make_topoO(
    // -------- 1-minute resolution
    blitz::Array<int16_t,2> const &FGICE1m,
    blitz::Array<int16_t,2> const &ZICETOP1m,
    blitz::Array<int16_t,2> const &ZSOLG1m,
    blitz::Array<int16_t,2> const &FOCEAN1m,
    // -------- 10-minute resolution
    blitz::Array<double,2> const &FLAKES)

{
    // ----------------------- Set up output variables
    ibmisc::ArrayBundle<double,2> out;
    auto &FOCEAN(out.add("FOCEAN", {
        "description", "0 or 1, Bering Strait 1 cell wide",
        "units", "1",
        "source", "GISS 1Qx1",
    }));
    auto &FLAKE(out.add("FLAKE", {
        "description", "Lake Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &FGRND(out.add("FGRND", {
        "description", "Ground Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &FGICE(out.add("FGICE", {
        "description", "Glacial Ice Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &FGICE_greenland(out.add("FGICE_greenland", {
        "description", "Greenland Ice Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &ZATMO(out.add("ZATMO", {
        "description", "Atmospheric Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &dZOCEN(out.add("dZOCEN", {
        "description", "Ocean Thickness",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &dZLAKE(out.add("dZLAKE", {
        "description", "Lake Thickness",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &dZGICE(out.add("dZGICE", {
        "description", "Glacial Ice Thickness",
        "units", "m",
        "sources", "Ekholm,Bamber",
    }));
    auto &ZSOLDG(out.add("ZSOLDG", {
        "description", "Solid Ground Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &ZICETOP(out.add("ZICETOP", {
        "description", "Atmospheric Topography (Ice-Covered Regions Only)",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &ZSGLO(out.add("ZSGLO", {
        "description", "Lowest Solid Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &ZLAKE(out.add("ZLAKE", {
        "description", "Lake Surface Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &ZGRND(out.add("ZGRND", {
        "description", "Topography Break between Ground and GIce",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &ZSGHI(out.add("ZSGHI", {
        "description", "Highest Solid Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    // ------------- Non-rounded versions
    auto &FOCEANF(out.add("FOCEANF", {
        "description", "Fractional ocean ocver",
        "units", "1",
        "sources", "GISS 1Qx1",
    }));
    auto &FGICEF(out.add("FGICEF", {
        "description", "Glacial Ice Surface Fraction (Ocean NOT rounded)",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &ZATMOF(out.add("ZATMOF", {
        "description", "Atmospheric Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));

    out.allocate({IM,JM}, {"im", "jm"},
        true, blitz::fortranArray);

    // ------------------------------------------------------
    double const TWOPI = 2. * M_PI;
    double const AREAG = 4. * M_PI;

    //
    // FOCEAN: Ocean Surface Fraction (0:1)
    //
    // Fractional ocean cover FOCEANF is interpolated from FOAAH2
    blitz::Array<double, 2> WT1m(const_array(shape(IM1m, JM1m), 1.0, FortranArray<2>()));
    Hntr hntr1q1m(17.17, g1qx1, g1mx1m);
    hntr1q1m.regrid<double,int16_t,double>(WT1m, FOCEAN1m, FOCEANF, true);    // Fractional ocean cover

    // --------- FGICE is interpolated from FGICE1m

    //hntr1q1m.regrid(FOCEAN1m, FGICE1m, FGICE, true, -1.0, 1.0);    // Use FCONT1m = 1-FOCEAN1m for weight
    // (use WT1m here for weight instead of 1-FOCEAN1m so that FOCEAN+FLAKE+FGICE = 1
    // (instead of FOCEAN+FLAKE+FGRND=1 and FGICE is a portion of FGRND).
    hntr1q1m.regrid(WT1m, FGICE1m, FGICEF, true);

    // Here, FGRND=1-FOCEAN is implied.

    int8_t const TO_NONE = 0;
    int8_t const TO_OCEAN = 1;
    int8_t const TO_LAND = 2;

    blitz::Array<int8_t,2> SET_TO(IM,JM, fortranArray);
    SET_TO = TO_NONE;

    // Set following grid cells to be continent

    // Antarctic Peninsula.  Gary Russell said (2018-02-20): "That
    // hole should be filled with land ice.  It is unclear how it got
    // there in the first place.  The ETOPO1 ocean fraction is bases
    // on ZICETOP being below 0 and adjacent to other ocean cells."
    SET_TO( 84, 18) = TO_LAND;
    SET_TO( 85, 18) = TO_LAND;
    SET_TO( 86, 18) = TO_LAND;

    SET_TO(236, 82) = TO_LAND;
    SET_TO(242, 82) = TO_LAND;
    SET_TO(245, 82) = TO_LAND;
    SET_TO(224,101) = TO_LAND;
    SET_TO( 53,119) = TO_LAND;
    SET_TO(171,125) = TO_LAND;  //  Cyprus
    SET_TO(164,126) = TO_LAND;  //  Crete
    SET_TO(158,129) = TO_LAND;
    SET_TO(158,130) = TO_LAND;
    SET_TO(242,131) = TO_LAND;
    SET_TO(263,136) = TO_LAND;
    SET_TO(258,137) = TO_LAND;
    SET_TO(258,138) = TO_LAND;
    SET_TO( 46,139) = TO_LAND;
    SET_TO(258,139) = TO_LAND;
    SET_TO(275,152) = TO_LAND;
    SET_TO(  8,156) = TO_LAND;
    SET_TO( 10,156) = TO_LAND;
    SET_TO( 12,157) = TO_LAND;
    SET_TO(172,157) = TO_LAND;
    SET_TO(202,157) = TO_LAND;
    SET_TO( 69,159) = TO_LAND;
    SET_TO(204,159) = TO_LAND;
    SET_TO( 62,167) = TO_LAND;
    SET_TO( 73,171) = TO_LAND;
    SET_TO( 75,171) = TO_LAND;
    SET_TO( 78,171) = TO_LAND;

    // Set following grid cells to be ocean
    SET_TO(179,105) = TO_OCEAN;
    SET_TO( 54,119) = TO_OCEAN;
    SET_TO(241,131) = TO_OCEAN;
    SET_TO(258,143) = TO_OCEAN;
    SET_TO(165,150) = TO_OCEAN;
    SET_TO(274,152) = TO_OCEAN;
    SET_TO( 15,154) = TO_OCEAN;
    SET_TO( 92,155) = TO_OCEAN;
    SET_TO( 13,157) = TO_OCEAN;
    SET_TO(173,157) = TO_OCEAN;
    SET_TO(176,157) = TO_OCEAN;
    SET_TO(203,157) = TO_OCEAN;
    SET_TO( 55,159) = TO_OCEAN;
    SET_TO(103,159) = TO_OCEAN;
    SET_TO(203,159) = TO_OCEAN;
    SET_TO( 67,160) = TO_OCEAN;
    SET_TO( 68,160) = TO_OCEAN;
    SET_TO( 79,160) = TO_OCEAN;
    SET_TO(199,160) = TO_OCEAN;
    SET_TO(126,161) = TO_OCEAN;
    SET_TO( 68,162) = TO_OCEAN;
    SET_TO( 75,165) = TO_OCEAN;
    SET_TO(225,169) = TO_OCEAN;  
    SET_TO(179,105) = TO_OCEAN;
    SET_TO( 54,119) = TO_OCEAN;
    SET_TO(241,131) = TO_OCEAN;
    SET_TO(258,143) = TO_OCEAN;
    SET_TO(165,150) = TO_OCEAN;
    SET_TO(274,152) = TO_OCEAN;
    SET_TO( 15,154) = TO_OCEAN;
    SET_TO( 92,155) = TO_OCEAN;
    SET_TO( 13,157) = TO_OCEAN;
    SET_TO(173,157) = TO_OCEAN;
    SET_TO(176,157) = TO_OCEAN;
    SET_TO(203,157) = TO_OCEAN;
    SET_TO( 55,159) = TO_OCEAN;
    SET_TO(103,159) = TO_OCEAN;
    SET_TO(203,159) = TO_OCEAN;
    SET_TO( 67,160) = TO_OCEAN;
    SET_TO( 68,160) = TO_OCEAN;
    SET_TO( 79,160) = TO_OCEAN;
    SET_TO(199,160) = TO_OCEAN;
    SET_TO(126,161) = TO_OCEAN;
    SET_TO( 68,162) = TO_OCEAN;
    SET_TO( 75,165) = TO_OCEAN;
    SET_TO(225,169) = TO_OCEAN;


    // South pole should be all ice; avoid singularities in regridding algo.
    for (int i=1; i<=IM; ++i) {
        FOCEAN(i,1) = 0;
        FGICE(i,1) = 1;
        FGICEF(i,1) = 1;
    }


    // FOCEAN (0 or 1) is rounded from FOCEANF
    for (int j=2; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        auto const foceanf(FOCEANF(i,j));
        bool to_ocean;
        switch(SET_TO(i,j)) {
            case TO_OCEAN :
                to_ocean = true;
                break;
            case TO_LAND :
                to_ocean = false;
                break;
            default:
                to_ocean = (foceanf >= .5);
        }

        if (to_ocean) {
            FOCEAN(i,j) = 1.;
            FGICE(i,j) = 0;
        } else {
            double const fact = 1. / (1. - foceanf);
//if (fact != 1) fprintf(stderr, "FACT (%d, %d): %g %g sum=%g, fact=%g\n", i,j, FGICE(i,j),  FOCEANF(i,j), FGICE(i,j)+FOCEANF(i,j), fact);
            FOCEAN(i,j) = 0.;
            FGICE(i,j) = FGICEF(i,j) * fact;
        }
    }}

#if 0
    // Average non-fractional and fractional ocean covers over latitude
    printf(
        " Comparison between Fractional and Non-fractional Ocean Cover\n\n"
        "         # of      # of     differ\n"
        "         fract    NOfrac      in #\n"
        "   J     cells     cells     cells\n"
        "   =     =====     =====     =====\n");
    blitz::Array<double,1> FOFLAT(JM, fortranArray);
    blitz::Array<double,1> FONLAT(JM, fortranArray);
    for (int J=JM; J >= 1; --J) {
        FOFLAT(J) = blitz::sum(FOCEANF(Range::all(),J));
        FONLAT(J) = blitz::sum(FOCEAN(Range::all(),J));
        printf("%4d%10.2f%10.2f%10.2f\n", J,FOFLAT(J),FONLAT(J),FOFLAT(J)-FONLAT(J));
    }
    double const factor = IM*JM / AREAG;

    double FOFSH = 0;
    double FONSH = 0;
    HntrGrid grid_g1qx1(g1qx1);
    for (int J=1; J <= JM/2; ++J) {
        FOFSH += factor * FOFLAT(J) * grid_g1qx1.dxyp(J);
        FONSH += factor * FONLAT(J) * grid_g1qx1.dxyp(J);
    }

    double FOFNH = 0;
    double FONNH = 0;
    for (int J=JM/2+1; J <= JM; ++J) {
        FOFNH += factor * FOFLAT(J) * grid_g1qx1.dxyp(J);
        FONNH += factor * FONLAT(J) * grid_g1qx1.dxyp(J);
    }

    printf("NH: %f %f %f\n", FOFNH, FONNH, FOFNH-FONNH);
    printf("SH: %f %f %f\n", FOFSH, FONSH, FOFSH-FONSH);
#endif

    //
    // FLAKE: Lake Surface Fraction (0:1)
    //
    // FLAKE is interpolated from FLAKES
    blitz::Array<double, 2> WTS(const_array(shape(IMS, JMS), 1.0, FortranArray<2>()));
    Hntr hntr1q10m(17.17, g1qx1, g10mx10m);
    hntr1q10m.regrid(WTS, FLAKES, FLAKE, true);

    // Antarctica and Arctic area have no lakes
    FLAKE(Range::all(), Range(1,JM/6)) = 0;             //  90:60 S
    FLAKE(Range::all(), Range(JM*14/15+1,JM)) = 0;      //  78:90 N
    FLAKE(Range(1,IM/2), Range(JM*41/45+1,JM)) = 0;  //  74:90 N, 0:180 W

    // TODO: Consider a more mask-based way to identify southern Greenland
    for (int J=(JM*5)/6; J <= (JM*11)/12; ++J) {    //  southern
    for (int I=IM/3+1; I <= (int)(.5 + .75*IM*(J-JM*.3)/JM); ++I) {  //  Greenland
       FLAKE(I,J) = 0;
    }}

    // Apportion FLAKE to the nonocean fraction and round to 1/256
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        double const fcont = 1. - FOCEAN(i,j);
        double const fcontf = 1. - FOCEANF(i,j);

        FLAKE(i,j) = FLAKE(i,j)*fcont / (fcontf+1e-20);
        // FLAKE(i,j) = round_mantissa_to(FLAKE(i,j), 8);
        FLAKE(i,j) = std::round(FLAKE(i,j)*256.) / 256.;
    }}

    //
    // FGICE: Glacial Ice Surface Fraction (0:1)
    //

    // Check that FGICE is between 0 and 1
    // If FGICE+FLAKE exceeds 1, reduce FLAKE
    for (int J=JM/6+1; J <= JM; ++J) {
    for (int I=1; I<=IM; ++I) {
        if (FGICE(I,J) < 0) {
            fprintf(stderr, "210: FGICE(%d,%d) < 0: %g\n" ,I,J,FGICE(I,J));
            FGICE(I,J) = 0;
        }
        if (FGICE(I,J) > 1) {
            fprintf(stderr, "210: FGICE(%d,%d) > 1: %g\n" ,I,J,FGICE(I,J));
            FGICE(I,J) = 1;
        }
        if (FLAKE(I,J)+FGICE(I,J)+FOCEAN(I,J) > 1) {
            fprintf(stderr, "210: FGICE+FLAKE+FOCEAN (%d,%d) > 1: %g + %g + %g\n",
                I,J,FGICE(I,J), FLAKE(I,J),FOCEAN(I,J));
            FLAKE(I,J) = 1. - FGICE(I,J) - FOCEAN(I,J);
        }
    }}

    // Replace land cells without vegetation with glacial ice in Z1QX1N
    // FGICE(35,52) = 1 - FLAKE(35,52)     4x3 Model

    //
    // FGRND: Surface Fraction of Ground (0:1)
    //

    // Check that FGRND is between 0 and 1
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        FGRND(i,j) = 1.0 - FOCEAN(i,j) - FLAKE(i,j) - FGICE(i,j);
        if (std::abs(FGRND(i,j)) < 1.e-14) FGRND(i,j) = 0;
        if (FGRND(i,j) < 0 || FGRND(i,j) > 1) {
            fprintf(stderr, "Error: FGRND(%d,%d) = %g %g %g %g\n", i,j,
                FGRND(i,j), FOCEAN(i,j),  FLAKE(i,j), FGICE(i,j));
        }
    }}

    //
    // dZOCEN: Ocean Thickness (m)
    //
    hntr1q1m.regrid(FOCEAN1m, ZSOLG1m, dZOCEN, true);
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        dZOCEN(i,j) = -dZOCEN(i,j) * FOCEAN(i,j);

        // Check that dZOCEN is positive
        if (FOCEAN(i,j) == 1 && dZOCEN(i,j) <= 0) {
            printf("Error: dZOCEN(%d,%d) <= 0 %g\n",i,j,dZOCEN(i,j));
        }
    }}

    //
    // dZGICE: Glacial Ice Thickness (m)
    //
    blitz::Array<double, 2> zictop(IM,JM, blitz::fortranArray);
    hntr1q1m.regrid(FGICE1m, ZICETOP1m, zictop, true);
    blitz::Array<double, 2> zsolg(IM,JM, blitz::fortranArray);
    hntr1q1m.regrid(FGICE1m, ZSOLG1m, zsolg, true);


    // RGICE = areal ratio of glacial ice to continent
    // For smaller ice caps and glaciers, dZGICH = CONSTK * RGICE^.3
    // Constant is chosen so that average value of dZGICH is 264.7 m
    // 264.7  =  sum(DXYP*FGIC1H*dZGICH) / sum(DXYP*FGIC1H)  =
    //        =  CONSTK * sum(DXYP*FGIC1H*RGICE^.3) / sum(DXYP*FGIC1H)
    dZGICE = 0;
    HntrGrid grid_g1qx1(g1qx1);
    {
        double SUMDFR = 0;
        double SUMDF = 0;
        blitz::Array<double,2> RGICE(IM, JM, blitz::fortranArray);
        for (int J=1; J<=JM; ++J) {
            double sum_dfr = 0;
            double sum_df = 0;
            for (int I=1; I<=IM; ++I) {
                double dzgice = zictop(I,J) - zsolg(I,J);

                if (dzgice != 0) {
                    // If this criterion (dzgice!=0) for GrIS+Ant doesn't work, use a mask instead.
                    // Set GrIS and AIS according to ETOPO1
                    dZGICE(I,J) = dzgice;
                } else {
                    double const fcont = 1. - FOCEAN(I,J);
                    RGICE(I,J) = FGICE(I,J) / (fcont+1e-20);
                    sum_dfr += FGICE(I,J) * std::pow(RGICE(I,J),.3);
                    sum_df += FGICE(I,J);
                }
            }
            SUMDFR += grid_g1qx1.dxyp(J) * sum_dfr;
            SUMDF += grid_g1qx1.dxyp(J) * sum_df;
        }

        double CONSTK = 264.7 * SUMDF / SUMDFR;

        // Replace in.FGICEH and dZGICH away from Greenland
        for (int J=1; J <= JM; ++J) {
        for (int I=1; I <= IM; ++I) {
            if (FGICE(I,J) == 0) {
                // Clear ice depth where there is no ice
                dZGICE(I,J) = 0;
            } else if (dZGICE(I,J) == 0) {
                // Avoid ice sheets
                dZGICE(I,J) = CONSTK * std::pow(RGICE(I,J), .3);
            }
        }}
    }

    //
    // ZATMO  = Atmospheric topography (m)
    // dZLAKE = Mean lake thickness (m)
    // ZSOLDG = Solid ground topography (m)
    // ZICETOP = Atmospheric topography, Ice covered regions (m)
    // ZSGLO  = Lowest value of ZICETOP1m in model cell (m)
    // ZLAKE  = Surface lake topography (m)
    // ZGRND  = Altitude break between ground and land ice (m)
    // ZSGHI  = Highest value of ZICETOP1m in model cell (m)
    //
    ZSOLDG = - dZOCEN;  //  solid ground topography of ocean
    ZICETOP = 0;
    callZ(
        FGICE1m, FOCEAN1m, ZICETOP1m, ZSOLG1m,
        FOCEAN, FLAKE, FGRND, ZATMO, ZATMOF,
        dZLAKE,ZSOLDG,ZICETOP,ZSGLO,ZLAKE,ZGRND,ZSGHI);

#if 0
This is ineffective: FOCEAN will be 0 or 1 already.

    // Adjust ZATMO so it's based on ENITRE gridcell,
    // not just land portion.
    for (int J=1; J <= JM; ++J) {
    for (int I=1; I<=IM; ++I) {
        ZATMO(I,J) = ZATMO(I,J) * (1 - FOCEAN(I,J));
//ZICETOP(I,J) -= ZATMO(I,J);
    }}
#endif


#if 1
    // Reset ZATMO, dZOCEN and ZLAKE by hand if FLAKE(I,J) == 1
    // This is here due to errors in ETOPO2
    // (Maybe) not needed for ETOPO1
    for (auto &reset : resets) {
        double elev = reset.elev;
        auto &ijs(reset.points);

        for (std::array<int,2> const &ij : ijs) {
            int const i(ij[0]);
            int const j(ij[1]);
            dZLAKE(i,j) += (elev - ZLAKE(i,j));
            ZATMO(i,j) = elev;
            ZLAKE(i,j) = elev;
        }
    }
#endif

    for (int J=1; J <= JM; ++J) {
    for (int I=1; I<=IM; ++I) {
        if (FLAKE(I,J) == 1) {
            fprintf(stderr, "FLAKE(%d,%d) == 1: %g %g %g\n",
                I, J, ZATMO(I,J),dZLAKE(I,J),ZLAKE(I,J));
        }
    }}


    // -------------------------------------------------
    // This will eventually be scaled to the Atmosphere (g2hx2) grid.
    // Adjust it on this grid so that the g2hx2 version is exactly
    // a regridded g1qx1 version.
    for (int jb=0; jb<JMB; ++jb) {
    for (int ib=0; ib<IMB; ++ib) {
        double foceanb = 0;
        for (int j=jb*2; j<jb*2+2; ++j) {
        for (int i=ib*2; i<ib*2+2; ++i) {
            foceanb += FOCEAN(i+1,j+1);
        }}

        if (foceanb > 0) {
            for (int j=jb*2; j<jb*2+2; ++j) {
                int const J=j+1;
                for (int i=ib*2; i<ib*2+2; ++i) {
                    int const I=i+1;
                    if (FLAKE(I,J) > 0) {
                        FGRND(I,J) += FLAKE(I,J);
                        FLAKE(I,J) = 0;
                        dZLAKE(I,J) = 0;
                    }
                }
            }
        }
    }}

    // Original code, which operate on the Atmosphere grid.
    // C**** Replace LAKES with GROUND if cell has some OCEAN
    // C****
    //       Do 10 J=1,JMB
    //       Do 10 I=1,IMB
    //       If (FOCEANB(I,J) > 0 .and. FLAKEB(I,J) > 0)  Then
    //            FGRNDB(I,J) = FGRNDB(I,J) + FLAKEB(I,J)
    //            FLAKEB(I,J) = 0
    //           dZLAKEB(I,J) = 0  ;  EndIf
    //    10 Continue

    return out;
}
// ---------------------------------------------------------------
// ======================================================================
MakeTopoO::MakeTopoO(
    FileLocator const &files,
    std::vector<std::string> const &_varinputs)
: hspec(*modele::grids.at("g1qx1"))
{
    // -------- 1-minute resolution
    blitz::Array<int16_t,2> FGICE1m(IM1m, JM1m, fortranArray);
    blitz::Array<int16_t,2> ZICETOP1m(IM1m, JM1m, fortranArray);
    blitz::Array<int16_t,2> ZSOLG1m(IM1m, JM1m, fortranArray);
    blitz::Array<int16_t,2> FOCEAN1m(IM1m, JM1m, fortranArray);
    // -------- 10-minute resolution (Z10MX10M.nc)
    blitz::Array<double,2> FLAKES(IMS, JMS, fortranArray);


    // Read into our variables, from user-specified locations
    NcBulkReader(&files, _varinputs)
        ("FGICE1m", FGICE1m)
        ("ZICETOP1m", ZICETOP1m)
        ("ZSOLG1m", ZSOLG1m)
        ("FOCEAN1m", FOCEAN1m)
        ("FLAKES", FLAKES);

printf("FINISHED READING INPUTS\n");

    bundle = _make_topoO(
        FGICE1m, ZICETOP1m, ZSOLG1m, FOCEAN1m,
        FLAKES);
}


}}
