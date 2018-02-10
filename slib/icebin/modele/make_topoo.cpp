#include <cmath>
#include <boost/filesystem.hpp>
#include <ibmisc/fortranio.hpp>
#include <icebin/error.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/make_topoo.hpp>


using namespace blitz;
using namespace ibmisc;
using namespace spsparse;

namespace icebin {
namespace modele {

static double const NaN = std::numeric_limits<double>::quiet_NaN();

static void callZ(
    // (IM1m, JM1m)
    blitz::Array<double,2> &FOCEAN1m,
    blitz::Array<double,2> &ZICETOP1m,
    blitz::Array<double,2> &ZSOLG1m,

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
    // Input:  FOCEAN1m = ocean fraction at 2 x 2 (minute)
    //         ZICETOP1m = solid topography (above ice) at 2 x 2 (minute)
    //         ZSOLG1m = solid ground topography at 2 x 2 (minute)
    //
    // Output: ZATMO  = atmospheric topography (m)
    //         dZLAKE = mean lake thickness (m)
    //         ZSOLDG = solid ground topography (m)
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
            int I1m = (IMAX == 1 ? IM1m : I*IM1m/IM);

            if (FOCEAN(I,J) != 0) {   // (I,J) is an ocean cell
                ZATMO(I,J) = 0;
                dZLAKE(I,J) = 0;
                // ZSOLDG(I,J) = - dZOCEN(I,J)  //  already filled in
                ZLAKE(I,J) = 0;
                ZGRND(I,J) = 0;
                ZSGHI(I,J) = 0;
                ZSGLO(I,J) = 999999;
                for (int J2=J11; J2 <= J1M; ++J2) {
                for (int I2=I11; I2 <= I1m; ++I2) {
                    if (ZSGLO(I,J) > ZICETOP1m(I2,J2) && FOCEAN1m(I2,J2) == 1.) {
                        ZSGLO(I,J) = ZICETOP1m(I2,J2);
                    }
                }}

                if (ZSGLO(I,J) == 999999) (*icebin_error)(-1,
                    "Ocean cell (%d,%d) has no ocean area on 2-minute grid; "
                    "I11,I1m,J11,J1M = %d,%d,%d,%d",
                    I,J,
                    I11,I1m,J11,J1M);
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
                for (int J1m=J11; J1m <= J1M; ++J1m) {
                for (int I2m=I11; I2m <= I1m; ++I2m) {
                    if (FOCEAN1m(I2m,J1m) == 1) continue;
                    double area = grid_g1mx1m.dxyp(J1m);
                    cells2.push_back(AreaDepth(area, ZICETOP1m(I2m,J1m)));

                    SAREA += area;
                    SAZSG += area*ZSOLG1m(I2m,J1m);
                }}
                std::sort(cells2.begin(), cells2.end());

                if (SAREA == 0) (*icebin_error)(-1,
                    "Continental cell (%d,%d) has no continental area on 2-minute grid. (%d-%d, %d-%d)",
                    I,J,I11,I1m,J11,J1M);

                // Determine ZSOLDG
                ZSOLDG(I,J) = SAZSG / SAREA;
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



blitz::ArrayBundle<double,2> make_topoO(
    // -------- 1-minute resolution
    blitz::Array<uint16_t,2> const &FGICE1m,
    blitz::Array<uint16_t,2> const &ZICETOP1m,
    blitz::Array<uint16_t,2> const &ZSOLG1m,
    blitz::Array<uint16_t,2> const &FOCEAN1m,
    // -------- 10-minute resolution
    blitz::Array<double,2> const &FLAKES)

{
    // ----------------------- Set up output variables
    blitz::ArrayBundle<double,2> out;
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
    auto &FOCENF(out.add("FOCENF", {
        "description", "Fractional ocean ocver",
        "units", "1",
        "sources", "GISS 1Qx1",
    }));
    out.allocate(blitz::shape(IM,JM), {"im", "jm"},
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
    hntr1q1m.regrid(WT1m, FOCEAN1m, FOCEANF, true);    // Fractional ocean cover

    // FOCEAN (0 or 1) is rounded from FOCEAN
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        FOCEAN(i,j) = std::round(FOCEANF(i,j));
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
        FLAKE(i,j) = FLAKE(i,j)*(1.-FOCEAN(i,j)) / (1.-FOCEANF(i,j)+1e-20);
        // FLAKE(i,j) = round_mantissa_to(FLAKE(i,j), 8);
        FLAKE(i,j) = std::round(FLAKE(i,j)*256.) / 256.;
    }}

    //
    // FGICE: Glacial Ice Surface Fraction (0:1)
    //
    // FGICE is interpolated from FGICE1m
    hntr1q1m.regrid(FCONT2, FGICE1m, FGICE, true);

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
    FGRND = 1 - FOCEAN - FLAKE - FGICE;

    // Check that FGRND is between 0 and 1
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        if (FGRND(i,j) < 0 || FGRND(i,j) > 1) {
            printf("Error: FGRND(%d,%d) = %g %g %g %g\n", i,j,
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
    hntr1q1m.regrid_cast<int16_t,double>(FGICE1m, ZICETOP1m, zictop, true);
    blitz::Array<double, 2> zsolg(IM,JM, blitz::fortranArray);
    hntr1q1m.regrid_cast<int16_t,double>(FGICE1m, ZSOLG1m, zsolg, true);

    for (int J=1; J<=JM; ++J) {
    for (int I=1; I<=IM; ++I) {
        double dzgice = zictop(I,J) - zsolg(I,J);
        if (FGICE(I,J) > 0) {
            dZGICE(I,J) = std::max(dzgice, 1.);
        } else {
            dZGICE(I,J) = 0;
        }
    }}

    //
    // ZATMO  = Atmospheric topography (m)
    // dZLAKE = Mean lake thickness (m)
    // ZSOLDG = Solid ground topography (m)
    // ZSGLO  = Lowest value of ZICETOP1m in model cell (m)
    // ZLAKE  = Surface lake topography (m)
    // ZGRND  = Altitude break between ground and land ice (m)
    // ZSGHI  = Highest value of ZICETOP1m in model cell (m)
    //
    ZSOLDG = - dZOCEN;  //  solid ground topography of ocean
    callZ(
        FOCEAN1m, ZICETOP1m, ZSOLG1m,
        FOCEAN, FLAKE, FGRND, ZATMO,
        dZLAKE,ZSOLDG,ZSGLO,ZLAKE,ZGRND,ZSGHI);

    // Reset ZATMO, dZOCEN and ZLAKE by hand if FLAKE(I,J) == 1
    //

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

    for (int J=1; J <= JM; ++J) {
    for (int I=1; I<=IM; ++I) {
        if (FLAKE(I,J) == 1) {
            printf("FLAKE(%d,%d) == 1: %g %g %g\n",
                I, J, ZATMO(I,J),dZLAKE(I,J),ZLAKE(I,J));
        }
    }}

    return out;
}
// ---------------------------------------------------------------
// ======================================================================
/** Reads a bunch of blitz::Arrays from a bunch of NetCDF files */
class BulkNcReader
{
    FileLocator const &files,
    std::map<std::string, std::array<std::string,2>> varmap;

    // Actions, keyed by filename
    std::vector<std::tuple<
        std::string, std::function<void (NcIO &ncio)>
    >> actions;

public:

    /** @param _vars {varname=interal name, fname=NetCDF filename, vname = NetcDF variable name, ...} */
    BulkNcReader(
        FileLocator const &_files,
        std::vector<std::string> const &_vars)
    : files(_files)
    {
        size_t n = _vars.size();
        if (3*(n/3) != n) (*icebin_error(-1,
            "BulkNcReader initializers must be in triplets: <varname>, <fname>, <vname>"));

        for (auto ii = _vars.begin(); ii != _vars.end(); ) {
            std::string const &var(*ii++);
            std::string const &fname(*ii++);
            std::string const &vname(*ii++);

            vars.insert(std::make_pair(var, make_array(fname, vname)));
        }
    }

private:

    template<class TypeT, int RANK>
    void add_var(
        blitz::Array<TypeT, RANK> &var,
        std::string const &fname,
        std::string const &vname)
    {
        actions.push_back(std::make_tuple(
            fname,
            std::bind(&ncio_blitz<TypeT,RANK>, _1, var, vname,
                "", std::vector<netCDF::NcDim>{})));

    }


public:

    /** @param varname Internal name for this variable, assigned in constructor */
    template<class TypeT, int RANK>
    void add_var(
        std::string const &varname,
        blitz::Array<TypeT, RANK> &var)
    {
        auto iix(varmap.find(varname));
        if (iix == varmap.end()) (*icebin_error)(-1,
            "User failed to include fname and vname for variable %s", varname.c_str());

        std::string const &fname(iix[0]);
        std::string connt &vname(iix[1]);

        add_var(var, fname, vname);
        varmap.erase(iix);    // We've added once, can't add again
    }


    void read_bulk()
    {
        // Check that every expected variable has been assigned.
        if (varmap.size() != 0) {
            for (auto ii=varmap.begin(); ii != varmap.end(); ++ii)
                fprintf("    Unassigned: %s\n", ii->first.c_str());
            (*icebin_error)(-1, "Unassigned variables");
        }


        // Read variables, grouped by file
        std::sort(actions.begin(), actions.end());
        std::unique_ptr<NcIO> ncio;
        std::string const last_fname = "";
        for (auto ii=actions.begin(); ii != actions.end(); ++ii) {
            auto &fname(std::get<0>(*ii));
            auto &action(std::get<1>(*ii));

            if (fname != last_fname) ncio.reset(new NcIO(files.locate(fname), 'r'));
            action(*ncio);
        }
    }

};


// ======================================================================
blitz::ArrayBundle<double,2> make_topoO(
    FileLocator const &files,
    std::vector<std::string> const &_varinputs)
{
    // -------- 1-minute resolution
    blitz::Array<uint16_t,2> FGICE1m(IM1m, JM1m, fortranArray),
    blitz::Array<uint16_t,2> ZICETOP1m(IM1m, JM1m, fortranArray),
    blitz::Array<uint16_t,2> ZSOLG1m(IM1m, JM1m, fortranArray),
    blitz::Array<uint16_t,2> FOCEAN1m(IM1m, JM1m, fortranArray),
    // -------- 10-minute resolution (Z10MX10M.nc)
    blitz::Array<double,2> FLAKES(IMS, JMS, fortranArray)


    // Read into our variables, from user-specified locations
    BulkNcReader bulk(_varinputs);
    bulk.add_var("FGICE1m", FGICE1m);
    bulk.add_var("ZICETOP1m", ZICETOP1m);
    bulk.add_var("ZSOLG1m", ZSOLG1m);
    bulk.add_var("FOCEAN1m", FOCEAN1m);
    bulk.add_var("FLAKES", FLAKES);
    bulk.read_bulk();

    return make_topoO(
        FGICE1m, ZICETOP1m, ZSOLG1m, FOCEAN1m,
        FLAKES);
}


}}
