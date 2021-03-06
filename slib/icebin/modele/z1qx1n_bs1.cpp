#include <cmath>
#include <boost/filesystem.hpp>
#include <ibmisc/fortranio.hpp>
#include <icebin/error.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/grids.hpp>
#include <icebin/modele/z1qx1n_bs1.hpp>

using namespace blitz;
using namespace ibmisc;
using namespace spsparse;

namespace icebin {
namespace modele {

static double const NaN = std::numeric_limits<double>::quiet_NaN();




// ==================================================================
ibmisc::ArrayBundle<int16_t,2> read_etopo1(std::string const &fname)
{
    ibmisc::ArrayBundle<int16_t,2> bundle;

    // --------- Arrays we return
    auto &focean(bundle.add("focean", {
        "units", "m",
        "description", "Ice Top Topography (dense indexing)"
    }));
    auto &zictop(bundle.add("zictop", {
        "units", "m",
        "description", "Ice Top Topography (dense indexing)"
    }));
    auto &zsolid(bundle.add("zsolid", {
        "units", "m",
        "description", "Solid Ground Topography (dense indexing)"
    }));
    bundle.allocate({IM1m, JM1m}, {"im1m", "jm1m"}, true, blitz::fortranArray);

    // ----------- Read into those arrays
    {fortran::UnformattedInput fin(fname, Endian::BIG);
    std::array<char,80> titlei;

        // Read from file
        fortran::read(fin) >> titlei >> focean >> fortran::endr;
        fortran::read(fin) >> titlei >> zictop >> fortran::endr;
        fortran::read(fin) >> titlei >> zsolid >> fortran::endr;
    }

    return bundle;
}

// ==================================================================
// -------------------------------------------------------
ArrayBundle<double,2> greenland_inputs_bundle(bool allocate)
{
    ArrayBundle<double, 2> bundle;
    bundle.add("FOCEN2", {IM2m, JM2m}, {"im2m", "jm2m"}, {});
    bundle.add("ZETOP2", {IM2m, JM2m}, {"im2m", "jm2m"}, {
        "description", "Solid Topography except for ice shelves",
        "units", "m",
        "source", "Z2MX2M.NGDC"
    });

    bundle.add("FOCENS", {IMS, JMS}, {"ims", "jms"}, {});

    bundle.add("FLAKES", {IMS, JMS}, {"ims", "jms"}, {});

    bundle.add("FGICEH", {IMH, JMH}, {"imh", "jmh"}, {});

    bundle.add("FCONT1", {IM1, JM1}, {"im1", "jm1"}, {});
    bundle.add("FGICE1", {IM1, JM1}, {"im1", "jm1"}, {});

    if (allocate) bundle.allocate(true, blitz::fortranArray);

    return bundle;
}

GreenlandInputs::GreenlandInputs(ArrayBundle<double,2> &&_bundle) :
    bundle(std::move(_bundle)),
    FOCEN2(bundle.array("FOCEN2")),
    ZETOP2(bundle.array("ZETOP2")),
    FOCENS(bundle.array("FOCENS")),
    FLAKES(bundle.array("FLAKES")),
    FGICEH(bundle.array("FGICEH")),
    FCONT1(bundle.array("FCONT1")),
    FGICE1(bundle.array("FGICE1"))
{
}



ArrayBundle<double,2> topo_inputs_bundle(bool allocate)
{

    ArrayBundle<double,2> bundle;
    bundle.add("FOCEN2", {IM2m, JM2m}, {"im2m", "jm2m"}, {
        "description", "Ocean Fraction",
        "units", "0 or 1",
    });
    bundle.add("ZETOP2", {IM2m, JM2m}, {"im2m", "jm2m"}, {
        "description", "Solid Topography except for ice shelves",
        "units", "m",
        "source", "Z2MX2M.NGDC"
    });

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
    FOCEN2(bundle.array("FOCEN2")),
    ZETOP2(bundle.array("ZETOP2")),
    FLAKES(bundle.array("FLAKES")),
    dZGICH(bundle.array("dZGICH")),
    FGICEH(bundle.array("FGICEH")),
    ZSOLDH(bundle.array("ZSOLDH")),
    FCONT1(bundle.array("FCONT1")),
    FGICE1(bundle.array("FGICE1"))
{
}




// ======================================================================
void read_raw(TopoInputs &in, bool separate, GreenlandInputs *greenland, FileLocator const &files)
{
    std::array<char,80> titlei;

    printf("BEGIN z1qx1n_bs1 Read Input Files\n");

    // ---------------------------------------------------------
    // Read in Z2MX2M.NGDC
    // Read hand-modified FOCEN2
    blitz::Array<double,2> greenland_focen2(
        greenland ? greenland->FOCEN2 : blitz::Array<double,2>(IM2m,JM2m,fortranArray));
    {NcIO ncio(files.locate("Z2MX2M.NGDC-SeparateGreenland.nc"), 'r');
        ncio_blitz(ncio, in.FOCEN2, "FOCEN2", "double",
            get_dims(ncio, {"jm2m", "im2m"}));    // ncdims in nc order

        ncio_blitz(ncio, in.ZETOP2, "ZETOP2", "double",
            get_dims(ncio, {"jm2m", "im2m"}));    // ncdims in nc order

        // Separate Greenland from the rest
        greenland_focen2 = NaN;

        for (int j=1; j<=JM2m; ++j) {
        for (int i=1; i<=IM2m; ++i) {
            if (in.FOCEN2(i,j) == 2.0) {
                if (separate) {
                    greenland_focen2(i,j) = 0.0;
                    if (greenland) greenland->ZETOP2(i,j) = in.ZETOP2(i,j);
                    in.FOCEN2(i,j) = 1.0;
                    in.ZETOP2(i,j) = -300.0;    // Greenland -> -300m ocean basin
                } else {
                    in.FOCEN2(i,j) = 0.0;
                }
            }
        }}
    }

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


    if (separate) {
        // Zero out lakes over Greenland
        blitz::Array<double, 2> WT2(const_array(shape(IM2m, JM2m), 1.0, FortranArray<2>()));
        Hntr hntr2mh(17.17, g10mx10m, g2mx2m, 0);
        blitz::Array<double,2> greenland_focens(
            greenland ? greenland->FOCENS : blitz::Array<double,2>(IMS,JMS,fortranArray));

        hntr2mh.regrid(WT2, greenland_focen2, greenland_focens, 0);

        if (greenland) greenland->FLAKES = NaN;
        for (int j=1; j<=JMS; ++j) {
        for (int i=1; i<=IMS; ++i) {
            if (std::abs(greenland_focens(i,j)) < 1e-14) {
                if (greenland) greenland->FLAKES(i,j) = in.FLAKES(i,j);
                in.FLAKES(i,j) = 0.0;
            }
        }}
    }

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
printf("BB1\n");
        ncio_blitz(ncio, in.FCONT1, "FCONT1", "double",
            get_dims(ncio, {"jm1", "im1"}));
printf("BB2\n");
        ncio_blitz(ncio, in.FGICE1, "FGICE1", "double",
            get_dims(ncio, {"jm1", "im1"}));
printf("BB3\n");

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


template<class TypeT, int RANK>
ibmisc::ArrayBundle<TypeT, RANK> elev_ice_bundle()
{
    ibmisc::ArrayBundle<TypeT, RANK> bundle;
    bundle.add("fgice", {
        "description", "Fractional ice cover (rounded to 0 or 1)",
        "units", "1",
        "source", "ETOPO1 + FGICEH"});
    bundle.add("elev", {
        "description", "Elevation, bottom of atmosphere",
        "units", "m",
        "source", "ETOPO1"});
    return bundle;
}

Etopo1Ice etopo1_ice(
    HntrSpec const &hntrA,              // Coarse grid
    blitz::Array<double,2> &FGICE,      // As read from Fortran file
    HntrSpec const &hntrI,              // Fine-scale grid
    blitz::Array<int16_t,2> &FOCEAN,      // As read from Fortran file
    blitz::Array<int16_t,2> &ZICTOP,    // As read from Fortran file
    blitz::Array<int16_t,2> &ZSOLID)
{
    Etopo1Ice ret(elev_ice_bundle<double,1>());
    auto &dimI(ret.dimI);

    blitz::Array<double,1> fgiceAs(reshape1(FGICE));
    blitz::Array<int16_t,1> foceanIs(reshape1(FOCEAN));
    blitz::Array<int16_t,1> zictopIs(reshape1(ZICTOP));
    blitz::Array<int16_t,1> zsolidIs(reshape1(ZSOLID));
    int nIs = foceanIs.extent(0);

    int const nIS = hntrI.size();



    // ---------------- Determine dimI

    // ------- Grid cells we know are 100% ice
    // Continental cells north of 78N are entirely glacial ice.
    for (int ijIs=(hntrI.jm*14/15) * hntrI.im; ijIs < nIs; ++ijIs)
        if (foceanIs(ijIs) == 0) dimI.add_dense(ijIs);

    // Add ETOP1 gridcells to the indexing
    // Count number of used grid cells
    for (int iIs=0; iIs<nIs; ++iIs)
        if (zictopIs(iIs) != zsolidIs(iIs)) dimI.add_dense(iIs);

    // Mark where the 100% ice grid cells end
    int const nId_Etopo1 = dimI.dense_extent();    // Indices that are JUST Etopo1

    // -------- Grid cells that may be <100% ice
    // Add to dimI indices from AvI (unscaled), derived from from fgiceA
    Hntr hntr(17.17, hntrA, hntrI);
    TupleList<int,double,2> AvI;

    struct IncludeIce {
        blitz::Array<double,1> const & fgice;    // 0-based
        IncludeIce(blitz::Array<double,1> const &_fgice) : fgice(_fgice) {}
        bool operator()(int i)  // true if we want to use this cell
            { return fgice(i) > 0.; }
    };
    hntr.overlap(
            accum::to_dense_ignore_missing(
                std::array<decltype(&dimI),2>{nullptr, &dimI},
            accum::ref(AvI)),
        1.0,    // Don't need radius here, we are scale free
        IncludeIce(fgiceAs));



    // ---------------------- Allocate
    int nId = dimI.dense_extent();

    blitz::Array<double,1> &fgiceId(ret.bundle.array("fgice"));
    blitz::Array<double,1> &elevId(ret.bundle.array("elev"));
    ret.bundle.allocate({nId}, {"nId"});
    fgiceId = NaN;
    elevId = NaN;

    // Densify into allocated space
    // Set 100% ice cells to 100% ice
    for (int iId=0; iId<nId_Etopo1; ++iId) {
        int const iIs = dimI.to_sparse(iId);
        fgiceId(iId) = 1.0;
        elevId(iId) = zictopIs(iIs);
    }
    // Fill in elevId for everywhere
    for (int iId = nId_Etopo1; iId < dimI.dense_extent(); ++iId) {
        int const iIs = dimI.to_sparse(iId);
        elevId(iId) = zictopIs(iIs);
    }

    // ---------- Determine snow-covered areas from fgiceA
    // If fgiceA(iA) = xyz%, then select the top xyz% of grid
    // cells in iI corresponding to iA.


    // sort AvI by A, then descending elevation in I...
    struct Compare {    // Tuple is (iAs, iId)
        blitz::Array<double,1> &elevId;
        Compare(blitz::Array<double,1> &_elevId) : elevId(_elevId) {}

        bool operator()(
            Tuple<int,double,2> const &a,
            Tuple<int,double,2> const &b)
        {
            if (a.index(0) < b.index(0)) return true;
            if (a.index(0) > b.index(0)) return false;
            return elevId(a.index(1)) > elevId(b.index(1));
        }
    };
    std::sort(AvI.tuples.begin(), AvI.tuples.end(), Compare(elevId));

    // Add a sentinel to make block-scan easier
    AvI.add({-1, -1}, 0.0);

    // Scan through by A grid cell (and the I cells that overlap it)
    HntrGrid hgridI(hntrI);    // Make dxyp available
    HntrGrid hgridA(hntrA);    // Make dxyp available
    int l0 = 0;
    int lastAs = AvI[0].index(0);
    for (int l1=1; l1<AvI.size(); ++l1) {
        auto const ijAs(AvI[l1].index(0));    // Sparse indexing
        auto const ijId(AvI[l1].index(1));    // Dense indexing

        if (ijAs != lastAs) {
            // ---------- Process block l0...l1

            // Get latitude index from 1-D index
            int jA = ijAs / hntrA.im;    // ModelE indexing

            // Snow-covered area for this A cell
            double snowA = fgiceAs(ijAs) * hgridA.dxyp(jA+1);     // [m^2]

            // Amount of snow-covered area we need to allocate to ice grid
            double remainI = snowA;
            for (int ix=l0; ix<l1; ++ix) {
                int const ijIs = dimI.to_sparse(ijId);    // to sparse indexing
                int const jI = ijIs / hntrI.im;    // ModelE 2-D indexing
                double const areaI = hgridI.dxyp(jI+1);

                if (areaI <= remainI) {
                    fgiceId(ijId) = 1.0;
                    remainI -= areaI;
                } else {
                    // I grid cell overflows allocated ice in A
                    fgiceId(ijId) = remainI / areaI;
                    remainI = 0;
                    break;
                }
            }
            lastAs = ijAs;
        }
    }

    // Round fgiceI to get final hi-res ice mask
    for (int iId=0; iId<fgiceId.extent(0); ++iId)
        fgiceId(iId) = std::round(fgiceId(iId));

    return ret;
}


void callZ(
    // (IM2m, JM2m)
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


{NcIO ncio("int.nc", 'w');

    ArrayBundle<double,2> bundle;

    // (IM2m, JM2m)
    bundle.add("FOCEN2", FOCEN2, {"im2m", "jm2m"}, {});
    bundle.add("ZSOLD2", ZSOLD2, {"im2m", "jm2m"}, {});
    bundle.add("ZSOLG2", ZSOLG2, {"im2m", "jm2m"}, {});

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

    HntrGrid grid_g2mx2m(g2mx2m);

    for (int J=1; J <= JM; ++J) {
        int J21 = (J-1)*JM2m/JM + 1;    // 2-minute cells inside (I,J)
        int J2M = J*JM2m/JM;
        int const IMAX= (J==1 || J==JM ? 1 : IM);
        for (int I=1; I<=IMAX; ++I) {
            int I21 = (I-1)*IM2m/IM + 1;
            int I2M = (IMAX == 1 ? IM2m : I*IM2m/IM);

            if (FOCEAN(I,J) != 0) {   // (I,J) is an ocean cell
                ZATMO(I,J) = 0;
                dZLAKE(I,J) = 0;
                // ZSOLDG(I,J) = - dZOCEN(I,J)  //  already filled in
                ZLAKE(I,J) = 0;
                ZGRND(I,J) = 0;
                ZSGHI(I,J) = 0;
                ZSGLO(I,J) = 999999;
                for (int J2=J21; J2 <= J2M; ++J2) {
                for (int I2=I21; I2 <= I2M; ++I2) {
                    if (ZSGLO(I,J) > ZSOLD2(I2,J2) && FOCEN2(I2,J2) == 1.) {
                        ZSGLO(I,J) = ZSOLD2(I2,J2);
                    }
                }}

                if (ZSGLO(I,J) == 999999) (*icebin_error)(-1,
                    "Ocean cell (%d,%d) has no ocean area on 2-minute grid; "
                    "I21,I2M,J21,J2M = %d,%d,%d,%d",
                    I,J,
                    I21,I2M,J21,J2M);
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
                for (int J2=J21; J2 <= J2M; ++J2) {
                for (int I2=I21; I2 <= I2M; ++I2) {
                    if (FOCEN2(I2,J2) == 1) continue;
                    double area = grid_g2mx2m.dxyp(J2);
                    cells2.push_back(AreaDepth(area, ZSOLD2(I2,J2)));

                    SAREA += area;
                    SAZSG += area*ZSOLG2(I2,J2);
                }}
                std::sort(cells2.begin(), cells2.end());

                if (SAREA == 0) (*icebin_error)(-1,
                    "Continental cell (%d,%d) has no continental area on 2-minute grid. (%d-%d, %d-%d)",
                    I,J,I21,I2M,J21,J2M);

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

    // ----------- Get hi-res ice mask and elevations
    ibmisc::ArrayBundle<int16_t,2> etopo1(read_etopo1(etopo1_fname));
    Etopo1Ice e1ice(etopo1_ice(
        ghxh,
            in.FGICEH,
        g1mx1m,
            etopo1.array("focean"),
            etopo1.array("zictop"),
            etopo1.array("zsolid")));
    etopo1.free();


#if 1
// ------- Test: convert what we have so far to 1-degree so we can see it.
{
    // ------------ Convert to 1-degree
    Hntr hntr(17.17, g1x1, g1mx1m);

    auto bundleA(elev_ice_bundle<double,2>());
    bundleA.allocate({g1x1.im, g1x1.jm}, {"im1", "jm1"}, true, blitz::fortranArray);

    struct MyAccum {
        SparseSet<int,int> &dimI;
        blitz::Array<double,2> &fgiceAs;
        blitz::Array<double,2> &elevAs;
        blitz::Array<double,1> &fgiceId;
        blitz::Array<double,1> &elevId;	// Mas is implied; masked-out cells aren't in dimI
        MyAccum(
            SparseSet<int,int> &_dimI,
            ArrayBundle<double,1> &bundleId,
            ArrayBundle<double,2> &bundleAs)
        :
            dimI(_dimI),
            fgiceId(bundleId.array("fgice")),
            elevId(bundleId.array("elev")),
            fgiceAs(bundleAs.array("fgice")),
            elevAs(bundleAs.array("elev"))
        {
            fgiceAs = 0;
            elevAs = 0;
        }
        void add(std::array<int,2> index, double val) {
            auto const iAs = index[0];
            auto const iIs = index[1];
            if (!dimI.in_sparse(iIs)) return;

            int const iId = dimI.to_dense(iIs);
            fgiceAs(iAs) += val * fgiceId(iId);
            elevAs(iAs) += val * elevId(iId);
        }
    };

    hntr.scaled_regrid_matrix(MyAccum(e1ice.dimI, e1ice.bundle, bundleA));


    // ------------ Store for manual checking
    NcIO ncio("allice.nc", 'w');
    bundleA.ncio(ncio, {}, "", "double");
}
return;
#endif


    //
    // Add small ice cap and glacier data to FGICEH and dZGICH
    // north of Antarctic area.

    blitz::Array<double, 2> WT1(const_array(blitz::shape(IM1, JM1), 1.0, FortranArray<2>()));
    Hntr hntr1h(17.17, ghxh, g1x1);
    auto FCON1H(hntr1h.regrid(WT1, in.FCONT1));
    auto FGIC1H(hntr1h.regrid(WT1, in.FGICE1));

    // RGIC1H = areal ratio of glacial ice to continent
    // For smaller ice caps and glaciers, dZGICH = CONSTK * RGIC1H^.3
    // Constant is chosen so that average value of dZGICH is 264.7 m
    // 264.7  =  sum(DXYP*FGIC1H*dZGICH) / sum(DXYP*FGIC1H)  =
    //        =  CONSTK * sum(DXYP*FGIC1H*RGIC1H^.3) / sum(DXYP*FGIC1H)
    double SUMDFR = 0;
    double SUMDF = 0;
    HntrGrid grid_ghxh(ghxh);
    blitz::Array<double,2> RGIC1H(IMH, JMH, fortranArray);
    for (int JH=1+JMH/6; JH<=JMH; ++JH) {
        double sum1 = 0;
        double sum2 = 0;
        for (int IH=1; IH<=IMH; ++IH) {
            if (in.FGICEH(IH,JH) > 0)
                FGIC1H(IH,JH) = 0;  //  ignore Greenland
            RGIC1H(IH,JH) = FGIC1H(IH,JH) / (FCON1H(IH,JH)+1e-20);
            sum1 += FGIC1H(IH,JH) * std::pow(RGIC1H(IH,JH),.3);
            sum2 += FGIC1H(IH,JH);
        }
        SUMDFR += grid_ghxh.dxyp(JH) * sum1;
        SUMDF += grid_ghxh.dxyp(JH) * sum2;
    }
    double CONSTK = 264.7 * SUMDF / SUMDFR;

    // Replace in.FGICEH and dZGICH away from Greenland
    for (int JH=JMH/6+1; JH <= JMH; ++JH) {
        for (int IH=1; IH <= IMH; ++IH) {
            if (in.FGICEH(IH,JH) == 0) {
                in.FGICEH(IH,JH) = FGIC1H(IH,JH);
                in.dZGICH(IH,JH) = CONSTK * std::pow(RGIC1H(IH,JH), .3);
            }
        }
    }


    // ETOPO2 treats Antarctic ice shelves as ocean.
    // When this happens ETOPO2 data are replaced with interpolated data
    // from in.FGICEH and dZGICH.  Resulting files are:
    // FOCEN2 = Ocean fraction (0 or 1) correct for Antarctic ice shelves
    // FCONT2 = Continental fraction (0 or 1)
    // FGICE2 = Glacial ice fraction (0 or 1)
    // dZGIC2 = Thickness of glacial ice (m)
    // ZSOLD2 = Solid topography (m)        (above ice)
    // ZSOLG2 = Solid ground topography (m) (beneath ice)
    //
    blitz::Array<double, 2> WTH(const_array(blitz::shape(IMH, JMH), 1.0, FortranArray<2>()));
    Hntr hntrhm2(17.17, g2mx2m, ghxh);
    auto FGICE2(hntrhm2.regrid(WTH, in.FGICEH));
    auto dZGIC2(hntrhm2.regrid(in.FGICEH, in.dZGICH));
    auto ZSOLD2(hntrhm2.regrid(in.FGICEH, in.ZSOLDH));

    // North of Antarctic area: 60S to 90N
    blitz::Array<double,2> FCONT2(IM2m, JM2m, fortranArray);
    blitz::Array<double,2> ZSOLG2(IM2m, JM2m, fortranArray);
    for (int J2=JM2m/6+1; J2 <= JM2m; ++J2) {
        FCONT2(Range::all(), J2) = 1. - in.FOCEN2(Range::all(), J2);
        FGICE2(Range::all(), J2) = FGICE2(Range::all(), J2) * FCONT2(Range::all(), J2);
        dZGIC2(Range::all(), J2) = dZGIC2(Range::all(), J2) * FCONT2(Range::all(), J2);
        ZSOLD2(Range::all(), J2) = in.ZETOP2(Range::all(), J2);
        ZSOLG2(Range::all(), J2) = in.ZETOP2(Range::all(), J2) - dZGIC2(Range::all(), J2);
    }

    // Antarctic area: 90S to 60S
    for (int J2=1; J2<=JM2m/6; ++J2) {
    for (int I2=1; I2<=IM2m; ++I2) {
        if (in.FOCEN2(I2,J2) == 0) {
            // in.FOCEN2(I2,J2) = 0
            FCONT2(I2,J2) = 1;
            FGICE2(I2,J2) = 1;
            // dZGIC2(I2,J2) = dZGIC2(I2,J2)  //  in.ZETOP2 has 2m and other low
            // ZSOLD2(I2,J2) = ZSOLD2(I2,J2)       values over ice shelves
            if (in.ZETOP2(I2,J2) >= 100)  ZSOLD2(I2,J2) = in.ZETOP2(I2,J2);
            ZSOLG2(I2,J2) = ZSOLD2(I2,J2) - dZGIC2(I2,J2);
        } else if (FGICE2(I2,J2) <= .5) {  //  and in.FOCEN2(I2,J2) == 1
            // in.FOCEN2(I2,J2) = 1
            FCONT2(I2,J2) = 0;
            FGICE2(I2,J2) = 0;
            dZGIC2(I2,J2) = 0;
            ZSOLD2(I2,J2) = in.ZETOP2(I2,J2);
            ZSOLG2(I2,J2) = in.ZETOP2(I2,J2);
        } else if (FGICE2(I2,J2) > .5) {  //  and in.FOCEN2(I2,J2) == 1
            in.FOCEN2(I2,J2) = 0;
            FCONT2(I2,J2) = 1;
            FGICE2(I2,J2) = 1;
            dZGIC2(I2,J2) = ZSOLD2(I2,J2) - in.ZETOP2(I2,J2);
            // ZSOLD2(I2,J2) = ZSOLD2(I2,J2)
            ZSOLG2(I2,J2) = in.ZETOP2(I2,J2);
        }
    }}


    //
    // FOCEAN: Ocean Surface Fraction (0:1)
    //
    // Fractional ocean cover FOCENF is interpolated from FOAAH2
    blitz::Array<double, 2> WT2(const_array(shape(IM2m, JM2m), 1.0, FortranArray<2>()));
    Hntr hntr2mq1(17.17, g1qx1, g2mx2m);
    hntr2mq1.regrid(WT2, in.FOCEN2, out.FOCENF, true);    // Fractional ocean cover

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
    // FGICE is interpolated from FGICE2
    hntr2mq1.regrid(FCONT2, FGICE2, out.FGICE, true);

    // Antarctica is entirely glacial ice, no lakes nor ground
    for (int j=1; j<=JM/6; ++j) {
        out.FGICE(Range::all(), j) = 1. - out.FOCEAN(Range::all(), j);
    }

    // Continental cells north of 78N are entirely glacial ice
    for (int j=(JM*14)/15+1; j<=JM; ++j) {
        out.FGICE(Range::all(), j) = 1. - out.FOCEAN(Range::all(), j);
    }

    // There is no glacial ice over oceans
    for (int j=JM/6+1; j<=JM; ++j) {
        out.FGICE(Range::all(), j) *= (1. - out.FOCEAN(Range::all(), j));
    }

    // Round FGICE to nearest 1/256
    for (int j=1; j<=JM; ++j) {
    for (int i=1; i<=IM; ++i) {
        out.FGICE(i,j) = std::round(out.FGICE(i,j)*256.) / 256.;
    }}



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
    hntr2mq1.regrid(in.FOCEN2, ZSOLG2, out.dZOCEN, true);
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
    hntr2mq1.regrid(FGICE2, dZGIC2, out.dZGICE, true);
    for (int J=1; J<=JM; ++J) {
    for (int I=1; I<=IM; ++I) {
        if (out.FGICE(I,J) > 0) {
            out.dZGICE(I,J) = std::max(out.dZGICE(I,J), 1.);
        } else {
            out.dZGICE(I,J) = 0;
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
    out.ZSOLDG = - out.dZOCEN;  //  solid ground topography of ocean
    callZ(
        in.FOCEN2,ZSOLD2,ZSOLG2,
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
