#include <iostream>
#include <everytrace.h>
#include <tclap/CmdLine.h>
#include <ibmisc/memory.hpp>
#include <ibmisc/fortranio.hpp>
#include <ibmisc/filesystem.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/grids.hpp>

using namespace ibmisc;
using namespace icebin::modele;


const int16_t GREENLAND_VAL = 2;    // Used to mark Greenland in ZNGDC1-SeparateGreenland/FCONT1

void store_h(NcIO &ncout, blitz::Array<int16_t,2> val1m, std::string const &vname)
{
    auto &valh(ncout.tmp.make<blitz::Array<double,2>>(JMH, IMH));
    auto valh_1(reshape1(valh));

    Hntr hntrh(17.17, ghxh, g1mx1m);

    auto val1m_1(reshape1(val1m));
    auto WTA(const_array(val1m_1.shape(), 1.0, blitz::FortranArray<1>()));
    hntrh.regrid_cast<double, int16_t, double>(WTA, val1m_1, valh_1, true);
    auto dimsh(get_or_add_dims(ncout, {"jmh", "imh"}, {JMH, IMH}));
    ncio_blitz(ncout, valh, vname, "double", dimsh);
}

// ============================================================
// ------ Files:
// * ETOPO1: Topo everywhere, Ice ONLY for GrIS and AIS
// * ZNGDC1-SeparateGreenland.nc: FCONT1=2 for Greenland.  Can be used
//     to remove GrIS from other datasets.  FGICE1 contains AIS, GrIS,
//     Alaska, Arctic, NOT Himalayas.
// * ZICEHXH: Only contains ice sheets.  Not useful for us
// * Z2MX2M.NGDC-SeparateGreenland.nc: FOCEN2=2 for Greenland.  Can be
//     used to remove GrIS from other datasets, at higher resolution

/** Produces a global map of ice extent and elevation on the etopo1 grid.
@param include_greenland Include Greenland ice in the ice map? */
void etopo1_ice(
    FileLocator const &files,
    bool include_greenland,
    std::string const &fname_out)
{
    // This is what we construct
    blitz::Array<int16_t,2> fgice1m(JM1m, IM1m);
    fgice1m = 0;

    // Read in ZNGDC1 (1-degree resolution)
    // Read hand-modified FCONT1 (formerly from "ZNGDC1")
    blitz::Array<double,2> fgice1(JM1,IM1);
    {NcIO ncio(files.locate("ZNGDC1-SeparateGreenland.nc"), 'r');
        // ------- Ice Fractions
        ncio_blitz(ncio, fgice1, "FGICE1", "double", {});
    }


    NcIO nch("globaliceh.nc", 'w');

    // ------- Process ETOPO1
    std::array<char,80> titlei;
    blitz::Array<int16_t,2> zicetop1m(JM1m, IM1m);
    blitz::Array<int16_t,2> focean1m(JM1m, IM1m);
    {std::string fname(files.locate("ZETOPO1.NCEI-SeparateGreenland.nc"));
        printf("Opening ETOPO1 at %s\n", fname.c_str());
        NcIO fin(files.locate("ZETOPO1.NCEI-SeparateGreenland.nc"));

        // Continental cells north of 78N are entirely glacial ice.
        // (Except for Greenland, if that is to be done on a separate ice sheet)
        {
            ncio_blitz(fin, focean1m, "FOCEAN", "short", {});

            // Continental cells north of 78N are entirely glacial ice.
            // (but ignore Greenland)
            for (int j1m=JM1m*14/15; j1m < JM1m; ++j1m) {
            for (int i1m=0; i1m < IM1m; ++i1m) {
                if (focean1m(j1m,i1m) == 0) {
                    // Convert indexing from 1minute --> 1degree grid
                    int const i1 = i1m / 60;
                    int const j1 = j1m / 60;

                    fgice1m(j1m, i1m) = 1;
                }
            }}
        }

        // Antarctica is supplied by ETOPO1
        {
            blitz::Array<int16_t,2> zsolid1m(JM1m, IM1m);

            // Read zicetop1m, zsolid1m
            ncio_blitz(fin, zicetop1m, "ZICTOP", "short", {});
            ncio_blitz(fin, zsolid1m, "ZSOLID", "short", {});

            // Use ETOPO1 for Southern Hemisphere Ice
            for (int j1m=0; j1m < JM1m/2; ++j1m) {
            for (int i1m=0; i1m < IM1m; ++i1m) {
                if ( (focean1m(j1m,i1m) == 0)
                    && (zicetop1m(j1m,i1m) != zsolid1m(j1m,i1m)))
                {
                    fgice1m(j1m, i1m) = 1;
                }
            }}

            // (Maybe) use EOTOPO1 for Greenland too
            if (include_greenland) {
                for (int j1m=JM1m/2; j1m < JM1m; ++j1m) {
                for (int i1m=0; i1m < IM1m; ++i1m) {
                    if ( (focean1m(j1m,i1m) == GREENLAND_VAL)
                        && (zicetop1m(j1m,i1m) != zsolid1m(j1m,i1m)))
                    {
                        fgice1m(j1m, i1m) = 1;
                    }
                }}
            }


            // Remove Greenland from zicetop1m, zsolid1m
            for (int j1m=JM1m/2; j1m < JM1m; ++j1m) {
            for (int i1m=0; i1m < IM1m; ++i1m) {
                if (focean1m(j1m,i1m) == GREENLAND_VAL) {
                    if (!include_greenland) {
                        zicetop1m(j1m,i1m) = -300;
                        zsolid1m(j1m,i1m) = -300;
                    }
                }}
            }

            // Store zicetop1m, zsolid1m
            {NcIO ncout(fname_out, 'w');
                auto dims(get_or_add_dims(ncout, {"jm1m", "im1m"}, {JM1m, IM1m}));

                ncio_blitz(ncout, zicetop1m, "ZICETOP1m", "short", dims);
                ncio_blitz(ncout, zsolid1m, "ZSOLG1m", "short", dims);

            }
            store_h(nch, zicetop1m, "ZICETOPh");
            store_h(nch, zsolid1m, "ZSOLGh");
        }
    }


    // Add northern-hemisphere non-Greenland ice specified in ZNGDC1
    // This must be downscaled to ETOPO1 grid
    auto g1mx1m_dxyp(make_dxyp(g1mx1m));
    auto g1x1_dxyp(make_dxyp(g1x1));
    std::vector<std::tuple<int16_t,int,int>> cells;
    for (int j1=JM1/2; j1 < JM1; ++j1) {
    for (int i1=0; i1 < IM1; ++i1) {
        double const snow1 = fgice1(j1, i1);
        if (snow1 == 0) continue;

        // Assemble cells in this gridcell, sorted by
        // descending elevation (increasing negative elevation)
        cells.clear();
        for (int j1m=j1*60; j1m < (j1+1)*60; ++j1m) {
        for (int i1m=i1*60; i1m < (i1+1)*60; ++i1m) {
            if (focean1m(j1m,i1m) == 0 || (include_greenland && focean1m(j1m,i1m)) == GREENLAND_VAL) {
                double const elev = zicetop1m(j1m,i1m);
                cells.push_back(std::make_tuple(-elev, j1m, i1m));
            }
        }}
        std::sort(cells.begin(), cells.end());
        if (cells.size() == 0) continue;

        // Snow-covered area for this cell in g1x1
        double remain1m = snow1 * g1x1_dxyp(j1);
printf("j1 i1=%d %d (area = %g %g)\n", j1, i1, snow1, remain1m);
        for (auto ii=cells.begin(); ii != cells.end(); ++ii) {
            int16_t const elev1m = -std::get<0>(*ii);
            int const j1m(std::get<1>(*ii));
            int const i1m(std::get<2>(*ii));
//printf("--> %d %d\n", j1m, i1m);
            double const area1m = g1mx1m_dxyp(j1m);
            if (area1m <= remain1m) {
                remain1m -= area1m;
                fgice1m(j1m, i1m) = 1;
            } else {
                // No landice on the sea

                // We're done; round the last grid cell on g1mx1m
                if (remain1m >= area1m * .5)
                    fgice1m(j1m, i1m) = 1;
                break;
            }
        }

    }}

    // Remove Greenland from focean1m
    for (int j1m=JM1m/2; j1m < JM1m; ++j1m) {
    for (int i1m=0; i1m < IM1m; ++i1m) {
        if (focean1m(j1m,i1m) == GREENLAND_VAL) {
            if (include_greenland) {
                focean1m(j1m,i1m) = 0;
            } else {
                focean1m(j1m,i1m) = 1;
            }
        }}
    }

    // Store zicetop1m, zsolid1m
    {NcIO ncout(fname_out, 'a');
        auto dims(get_or_add_dims(ncout, {"jm1m", "im1m"}, {JM1m, IM1m}));
        ncio_blitz(ncout, focean1m, "FOCEAN1m", "short", dims);
        ncio_blitz(ncout, fgice1m, "FGICE1m", "short", dims);

    }
    store_h(nch, focean1m, "FOCEANh");
    store_h(nch, fgice1m, "FGICEh");

}

// ============================================================

struct ParseArgs {
    std::string ofname;
    bool greenland;

    ParseArgs(int argc, char **argv);
};

ParseArgs::ParseArgs(int argc, char **argv)
{
    // Wrap everything in a try block.  Do this every time, 
    // because exceptions will be thrown for problems.
    try {  
        TmpAlloc tmp;

        // Define the command line object, and insert a message
        // that describes the program. The "Command description message" 
        // is printed last in the help text. The second argument is the 
        // delimiter (usually space) and the last one is the version number. 
        // The CmdLine object parses the argv array based on the Arg objects
        // that it contains. 
        TCLAP::CmdLine cmd("Command description message", ' ', "<no-version>");

        TCLAP::UnlabeledValueArg<std::string> ofname_a(
            "ofname", "Name of output file", true, "etopo1_ice.nc", "output filename", cmd);
        TCLAP::SwitchArg greenland_a("g", "greenland", "Include Greenland?", cmd, false);

        // Parse the argv array.
        cmd.parse( argc, argv );

        // Get the value parsed by each arg.
        ofname = ofname_a.getValue();
        greenland = greenland_a.getValue();

    } catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        exit(1);
    }
}


int main(int argc, char** argv)
{
    using namespace TCLAP;

    everytrace_init();
    ParseArgs args(argc, argv);

    // Read the input files
    etopo1_ice(EnvSearchPath("MODELE_FILE_PATH"), args.greenland, args.ofname);
}


#if 0
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

