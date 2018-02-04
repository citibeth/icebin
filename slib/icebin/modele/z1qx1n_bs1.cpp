#include <cmath>
#include <boost/filesystem.hpp>
#include <ibmisc/fortranio.hpp>
#include <icebin/error.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/z1qx1n_bs1.hpp>

using namespace blitz;
using namespace ibmisc;
using namespace spsparse;

namespace icebin {
namespace modele {

static double const NaN = std::numeric_limits<double>::quiet_NaN();



ArrayBundle<double,2> topo_inputs_bundle(bool allocate)
{

    ArrayBundle<double,2> bundle;
    // ---------------------------------------------------
    // ---- Output from etopo1_ice
    bundle.add("FOCEN1m", {IM1m, JM1m}, {"im2", "jm2"}, {
        "description", "Ocean Fraction",
    });
    bundle.add("FGICE1m", {IM1m, JM1m}, {"im2", "jm2"}, {
        "description", "Ocean Fraction",
    });

    bundle.add("ZICTOP1m", {IM1m, JM1m}, {"im2", "jm2"}, {
        "description", "Solid Topography, top of ice",
        "units", "m",
    });
    bundle.add("ZSOLG1m", {IM1m, JM1m}, {"im2", "jm2"}, {
        "description", "Solid Bedrock Topography",
        "units", "m",
    });

    // ---------------------------------------------------
    bundle.add("FLAKES", {IMS, JMS}, {"ims", "jms"}, {
        "description", "Lake Fraction",
        "units", "0:1",
        "source", "Z2MX2M.NGDC"
    });

    bundle.add("dZGICH", {IMH, JMH}, {"imh", "jmh"}, {
        "description", "Glacial Ice Thickness",
        "units", "m",
        "source", "ZICEHXH"
    });
    bundle.add("FGICEH", {IMH, JMH}, {"imh", "jmh"}, {
        "description", "Glacial Ice Fraction (Antarctica & Greenland only)",
        "units", "0:1",
        "source", "ZICEHXH"
    });
    bundle.add("ZSOLDH", {IMH, JMH}, {"imh", "jmh"}, {
        "description", "Ice Topography (Antarctica & Greenland only)",
        "units", "m",
        "source", "ZICEHXH"
    });

    bundle.add("FCONT1", {IM1, JM1}, {"im1", "jm1"}, {
        "description", "Continental Fraction",
        "units", "0:1",
        "SOURCE", "ZNGDC1"
    });
    bundle.add("FGICE1", {IM1, JM1}, {"im1", "jm1"}, {
        "description", "Glacial Ice Fraction (all ice)",
        "units", "0:1",
        "SOURCE", "ZNGDC1"
    });

    if (allocate) {
        bundle.allocate(true, blitz::fortranArray);
    }
    return bundle;
}

TopoInputs::TopoInputs(ArrayBundle<double,2> &&_bundle) :
    bundle(std::move(_bundle)),
    FOCEN1m(bundle.array("FOCEN1m")),
    FGICE1m(bundle.array("FGICE1m")),
    ZICTOP1m(bundle.array("ZICTOP1m")),
    ZSOLG1m(bundle.array("ZSOLG1m")),
   

    FLAKES(bundle.array("FLAKES")),
    dZGICH(bundle.array("dZGICH")),
    FGICEH(bundle.array("FGICEH")),
    ZSOLDH(bundle.array("ZSOLDH")),
    FCONT1(bundle.array("FCONT1")),
    FGICE1(bundle.array("FGICE1"))
{
}




// ======================================================================
void read_raw(TopoInputs &in, FileLocator const &files)
{
    std::array<char,80> titlei;

    printf("BEGIN z1qx1n_bs1 Read Input Files\n");

    // ---------------------------------------------------------
    // Read output of etopo1_ice (on 2mx2m grid)
    // Read hand-modified FOCEN1m
    TODO

    // ---------------------------------------------------------
    // Read in Z10MX10M
    {std::string const fname = "Z10MX10M";
        fortran::UnformattedInput fin(files.locate(fname), Endian::BIG);
        fortran::read(fin) >> fortran::endr;
        fortran::read(fin) >> titlei >>
            fortran::blitz_cast<float, double, 2>(in.FLAKES) >> fortran::endr;
        //in.bundle.at("FLAKES").description = fortran::trim(titlei);
        printf("FLAKES read from %s: %s\n", fname.c_str(), fortran::trim(titlei).c_str());
    }


    // Zero out lakes over Greenland
#if 0
TODO
    blitz::Array<double, 2> WT2(const_array(shape(IM1m, JM1m), 1.0, FortranArray<2>()));
    Hntr hntr1mh(17.17, g10mx10m, g1mx1m, 0);
    blitz::Array<double,2> greenland_focens(
        greenland ? greenland->FOCENS : blitz::Array<double,2>(IMS,JMS,fortranArray));

    hntr1mh.regrid(WT2, greenland_focen2, greenland_focens, 0);

    if (greenland) greenland->FLAKES = NaN;
    for (int j=1; j<=JMS; ++j) {
    for (int i=1; i<=IMS; ++i) {
        if (std::abs(greenland_focens(i,j)) < 1e-14) {
            if (greenland) greenland->FLAKES(i,j) = in.FLAKES(i,j);
            in.FLAKES(i,j) = 0.0;
        }
    }}
#endif

    // ---------------------------------------------------------
    // Read in ZICEHXH
    {std::string const fname = "ZICEHXH";
        fortran::UnformattedInput fin(files.locate(fname), Endian::BIG);

        fortran::read(fin) >> titlei >>
            fortran::blitz_cast<float, double, 2>(in.dZGICH) >> fortran::endr;
        //in.bundle.at("dZGICH").description = fortran::trim(titlei);
        printf("dZGICH read from %s: %s\n", fname.c_str(), fortran::trim(titlei).c_str());

        fortran::read(fin) >> titlei >>
            fortran::blitz_cast<float, double, 2>(in.FGICEH) >> fortran::endr;
        //in.bundle.at("FGICEH").description = fortran::trim(titlei);
        printf("FGICEH read from %s: %s\n", fname.c_str(), fortran::trim(titlei).c_str());

        fortran::read(fin) >> titlei >>
            fortran::blitz_cast<float, double, 2>(in.ZSOLDH) >> fortran::endr;
        //in.bundle.at("ZSOLDH").description = fortran::trim(titlei);
        printf("ZSOLDH read from %s: %s\n", fname.c_str(), fortran::trim(titlei).c_str());
    }

    // Separate Greenland
    if (separate) {
        if (greenland) greenland->FGICEH = NaN;
        for (int j=JMH*3/4+1; j<=JMH; ++j) {
        for (int i=1; i<=IMH; ++i) {
            if (in.FGICEH(i,j) != 0) {
                if (greenland) greenland->FGICEH(i,j) = in.FGICEH(i,j);
                in.FGICEH(i,j) = 0;
            }
        }}
    }

    // -------------------------------------------------------------------
    // Read in ZNGDC1
    // Read hand-modified FCONT1 (formerly from "ZNGDC1")
    {NcIO ncio(files.locate("ZNGDC1-SeparateGreenland.nc"), 'r');
        ncio_blitz(ncio, in.FCONT1, "FCONT1", "double",
            get_dims(ncio, {"jm1", "im1"}));
        ncio_blitz(ncio, in.FGICE1, "FGICE1", "double",
            get_dims(ncio, {"jm1", "im1"}));

        // Separate Greenland from the rest
        if (greenland) {
            greenland->FCONT1 = NaN;
            greenland->FGICE1 = NaN;
        }
        for (int j=1; j<=JM1; ++j) {
        for (int i=1; i<=IM1; ++i) {
            if (in.FCONT1(i,j) == 2.0) {
                if (separate) {
                    if (greenland) {
                        greenland->FCONT1(i,j) = 1.0;
                        greenland->FGICE1(i,j) = in.FGICE1(i,j);
                    }
                    in.FCONT1(i,j) = 0;
                    in.FGICE1(i,j) = 0;
                } else {
                    in.FCONT1(i,j) = 1.0;
                }
            }
        }}
    }
    printf("END z1qx1n_bs1 Read Input Files\n");
}


void callZ(
    // (IM1m, JM1m)
    blitz::Array<double,2> &FOCEN1m,
    blitz::Array<double,2> &ZICTOP1m,
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


{NcIO ncio("int.nc", 'w');

    ArrayBundle<double,2> bundle;

    // (IM1m, JM1m)
    bundle.add("FOCEN1m", FOCEN1m, {"im2", "jm2"}, {});
    bundle.add("ZICTOP1m", ZICTOP1m, {"im2", "jm2"}, {});
    bundle.add("ZSOLG1m", ZSOLG1m, {"im2", "jm2"}, {});

    // (IM, IM)
    bundle.add("FOCEAN", FOCEAN, {"im", "jm"}, {});
    bundle.add("FLAKE", FLAKE, {"im", "jm"}, {});
    bundle.add("FGRND", FGRND, {"im", "jm"}, {});

    // (IM, IM)
    bundle.add("ZATMO", ZATMO, {"im", "jm"}, {});
    bundle.add("dZLAKE", dZLAKE, {"im", "jm"}, {});
    bundle.add("ZSOLDG", ZSOLDG, {"im", "jm"}, {});
    bundle.add("ZSGLO", ZSGLO, {"im", "jm"}, {});
    bundle.add("ZLAKE", ZLAKE, {"im", "jm"}, {});
    bundle.add("ZGRND", ZGRND, {"im", "jm"}, {});
    bundle.add("ZSGHI", ZSGHI, {"im", "jm"}, {});

    bundle.ncio(ncio, {}, "", "double");
}



    //
    // Input:  FOCEN1m = ocean fraction at 2 x 2 (minute)
    //         ZICTOP1m = solid topography (above ice) at 2 x 2 (minute)
    //         ZSOLG1m = solid ground topography at 2 x 2 (minute)
    //
    // Output: ZATMO  = atmospheric topography (m)
    //         dZLAKE = mean lake thickness (m)
    //         ZSOLDG = solid ground topography (m)
    //         ZSGLO  = lowest value of ZICTOP1m in model cell (m)
    //         ZLAKE  = surface lake topography (m)
    //         ZGRND  = altitude break between ground and land ice (m)
    //         ZSGHI  = highes value of ZICTOP1m in model cell (m)
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
                    if (ZSGLO(I,J) > ZICTOP1m(I2,J2) && FOCEN1m(I2,J2) == 1.) {
                        ZSGLO(I,J) = ZICTOP1m(I2,J2);
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
                    if (FOCEN1m(I2m,J1m) == 1) continue;
                    double area = grid_g1mx1m.dxyp(J1m);
                    cells2.push_back(AreaDepth(area, ZICTOP1m(I2m,J1m)));

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



void z1qx1n_bs1(TopoInputs &in, std::string const &etopo1_fname, TopoOutputs<2> &out)
{
    double const TWOPI = 2. * M_PI;
    double const AREAG = 4. * M_PI;


    //
    // FOCEAN: Ocean Surface Fraction (0:1)
    //
    // Fractional ocean cover FOCENF is interpolated from FOAAH2
    blitz::Array<double, 2> WT2(const_array(shape(IM1m, JM1m), 1.0, FortranArray<2>()));
    Hntr hntr1mq1(17.17, g1qx1, g1mx1m);
    hntr1mq1.regrid(WT2, in.FOCEN1m, out.FOCENF, true);    // Fractional ocean cover

    // FOCEAN (0 or 1) is rounded from FOCEAN
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        out.FOCEAN(i,j) = std::round(out.FOCENF(i,j));
    }}

    // Set following grid cells to be continent
    out.FOCEAN( 84, 18) = 0;
    out.FOCEAN( 85, 18) = 0;
    out.FOCEAN(236, 82) = 0;
    out.FOCEAN(242, 82) = 0;
    out.FOCEAN(245, 82) = 0;
    out.FOCEAN(224,101) = 0;
    out.FOCEAN( 53,119) = 0;
    out.FOCEAN(171,125) = 0;  //  Cyprus
    out.FOCEAN(164,126) = 0;  //  Crete
    out.FOCEAN(158,129) = 0;
    out.FOCEAN(158,130) = 0;
    out.FOCEAN(242,131) = 0;
    out.FOCEAN(263,136) = 0;
    out.FOCEAN(258,137) = 0;
    out.FOCEAN(258,138) = 0;
    out.FOCEAN( 46,139) = 0;
    out.FOCEAN(258,139) = 0;
    out.FOCEAN(275,152) = 0;
    out.FOCEAN(  8,156) = 0;
    out.FOCEAN( 10,156) = 0;
    out.FOCEAN( 12,157) = 0;
    out.FOCEAN(172,157) = 0;
    out.FOCEAN(202,157) = 0;
    out.FOCEAN( 69,159) = 0;
    out.FOCEAN(204,159) = 0;
    out.FOCEAN( 62,167) = 0;
    out.FOCEAN( 73,171) = 0;
    out.FOCEAN( 75,171) = 0;
    out.FOCEAN( 78,171) = 0;

    // Set following grid cells to be ocean
    out.FOCEAN(179,105) = 1;
    out.FOCEAN( 54,119) = 1;
    out.FOCEAN(241,131) = 1;
    out.FOCEAN(258,143) = 1;
    out.FOCEAN(165,150) = 1;
    out.FOCEAN(274,152) = 1;
    out.FOCEAN( 15,154) = 1;
    out.FOCEAN( 92,155) = 1;
    out.FOCEAN( 13,157) = 1;
    out.FOCEAN(173,157) = 1;
    out.FOCEAN(176,157) = 1;
    out.FOCEAN(203,157) = 1;
    out.FOCEAN( 55,159) = 1;
    out.FOCEAN(103,159) = 1;
    out.FOCEAN(203,159) = 1;
    out.FOCEAN( 67,160) = 1;
    out.FOCEAN( 68,160) = 1;
    out.FOCEAN( 79,160) = 1;
    out.FOCEAN(199,160) = 1;
    out.FOCEAN(126,161) = 1;
    out.FOCEAN( 68,162) = 1;
    out.FOCEAN( 75,165) = 1;
    out.FOCEAN(225,169) = 1;

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
        FOFLAT(J) = blitz::sum(out.FOCENF(Range::all(),J));
        FONLAT(J) = blitz::sum(out.FOCEAN(Range::all(),J));
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

    //
    // FLAKE: Lake Surface Fraction (0:1)
    //
    // FLAKE is interpolated from FLAKES
    blitz::Array<double, 2> WTS(const_array(shape(IMS, JMS), 1.0, FortranArray<2>()));
    Hntr hntr10m1q(17.17, g1qx1, g10mx10m);
    hntr10m1q.regrid(WTS, in.FLAKES, out.FLAKE, true);

    // Antarctica and Arctic area have no lakes
    out.FLAKE(Range::all(), Range(1,JM/6)) = 0;             //  90:60 S
    out.FLAKE(Range::all(), Range(JM*14/15+1,JM)) = 0;      //  78:90 N
    out.FLAKE(Range(1,IM/2), Range(JM*41/45+1,JM)) = 0;  //  74:90 N, 0:180 W

    for (int J=(JM*5)/6; J <= (JM*11)/12; ++J) {    //  southern
    for (int I=IM/3+1; I <= (int)(.5 + .75*IM*(J-JM*.3)/JM); ++I) {  //  Greenland
       out.FLAKE(I,J) = 0;
    }}

    // Apportion out.FLAKE to the nonocean fraction and round to 1/256
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        out.FLAKE(i,j) = out.FLAKE(i,j)*(1.-out.FOCEAN(i,j)) / (1.-out.FOCENF(i,j)+1e-20);
        out.FLAKE(i,j) = std::round(out.FLAKE(i,j)*256.) / 256.;
    }}

    //
    // FGICE: Glacial Ice Surface Fraction (0:1)
    //
    // FGICE is interpolated from FGICE1m
    hntr1mq1.regrid(FCONT2, in.FGICE1m, out.FGICE, true);

    // Check that FGICE is between 0 and 1
    // If out.FGICE+FLAKE exceeds 1, reduce FLAKE
    for (int J=JM/6+1; J <= JM; ++J) {
    for (int I=1; I<=IM; ++I) {
        if (out.FGICE(I,J) < 0) {
            fprintf(stderr, "210: out.FGICE(%d,%d) < 0: %g\n" ,I,J,out.FGICE(I,J));
            out.FGICE(I,J) = 0;
        }
        if (out.FGICE(I,J) > 1) {
            fprintf(stderr, "210: out.FGICE(%d,%d) > 1: %g\n" ,I,J,out.FGICE(I,J));
            out.FGICE(I,J) = 1;
        }
        if (out.FLAKE(I,J)+out.FGICE(I,J)+out.FOCEAN(I,J) > 1) {
            fprintf(stderr, "210: FGICE+FLAKE+out.FOCEAN (%d,%d) > 1: %g + %g + %g\n",
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
    hntr1mq1.regrid(in.FOCEN1m, in.ZSOLG1m, out.dZOCEN, true);
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        out.dZOCEN(i,j) = -out.dZOCEN(i,j) * out.FOCEAN(i,j);

        // Check that out.dZOCEN is positive
        if (out.FOCEAN(i,j) == 1 && out.dZOCEN(i,j) <= 0) {
            printf("Error: out.dZOCEN(%d,%d) <= 0 %g\n",i,j,out.dZOCEN(i,j));
        }
    }}

    //
    // dZGICE: Glacial Ice Thickness (m)
    //
    blitz::Array<double, 2> zictop(IM,JM, blitz::fortranArray);
    hntr1mq1.regrid_cast<int16_t,double>(in.FGICE1m, in.ZICTOP1m, zictop, true);
    blitz::Array<double, 2> zsolg(IM,JM, blitz::fortranArray);
    hntr1mq1.regrid_cast<int16_t,double>(in.FGICE1m, in.ZSOLG1m, zsolg, true);

    for (int J=1; J<=JM; ++J) {
    for (int I=1; I<=IM; ++I) {
        double dzgice = zictop(I,J) - zsolg(I,J);
        if (out.FGICE(I,J) > 0) {
            out.dZGICE(I,J) = std::max(dzgice, 1.);
        } else {
            out.dZGICE(I,J) = 0;
        }
    }}

    //
    // ZATMO  = Atmospheric topography (m)
    // dZLAKE = Mean lake thickness (m)
    // ZSOLDG = Solid ground topography (m)
    // ZSGLO  = Lowest value of ZICTOP1m in model cell (m)
    // ZLAKE  = Surface lake topography (m)
    // ZGRND  = Altitude break between ground and land ice (m)
    // ZSGHI  = Highest value of ZICTOP1m in model cell (m)
    //
    out.ZSOLDG = - out.dZOCEN;  //  solid ground topography of ocean
    callZ(
        in.FOCEN1m,in.ZICTOP1m,in.ZSOLG1m,
        out.FOCEAN, out.FLAKE, out.FGRND, out.ZATMO,
        out.dZLAKE,out.ZSOLDG,out.ZSGLO,out.ZLAKE,out.ZGRND,out.ZSGHI);

    // Reset ZATMO, out.dZOCEN and ZLAKE by hand if FLAKE(I,J) == 1
    //

    for (auto &reset : resets) {
        double elev = reset.elev;
        auto &ijs(reset.points);

        for (std::array<int,2> const &ij : ijs) {
            int const i(ij[0]);
            int const j(ij[1]);
            out.dZLAKE(i,j) += (elev - out.ZLAKE(i,j));
            out.ZATMO(i,j) = elev;
            out.ZLAKE(i,j) = elev;
        }
    }

    for (int J=1; J <= JM; ++J) {
    for (int I=1; I<=IM; ++I) {
        if (out.FLAKE(I,J) == 1) {
            printf("FLAKE(%d,%d) == 1: %g %g %g\n",
                I, J, out.ZATMO(I,J),out.dZLAKE(I,J),out.ZLAKE(I,J));
        }
    }}
}


}}
