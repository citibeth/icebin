#include <iostream>
#include <tclap/CmdLine.h>
#include <icebin/modele/z1qx1n_bs1.hpp>
#include <ibmisc/memory.hpp>
#include <everytrace.h>

using namespace ibmisc;
using namespace icebin::modele;



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
    bool include_greenland)
{
    // This is what we construct
    blitz::Array<int16_t,2> fgice1m(JM1m, IM1m);
    fgice1m = 0;

    // Read in ZNGDC1 (1-degree resolution)
    // Read hand-modified FCONT1 (formerly from "ZNGDC1")
    blitz::Array<int16_t,2> fcont1(JM1,IM1);    // cast to int16_t
    blitz::Array<double,2> fgice1(JM1,IM1);
    {NcIO ncio(files.locate("ZNGDC1-SeparateGreenland.nc"), 'r');
        // ------- GrIS Mask
        ncio_blitz(ncio, fcont1, "FCONT1", "double", {});
        // ------- Ice Fractions
        ncio_blitz(ncio, fgice1, "FGICE1", "double", {});
    }


    // ------- Process ETOPO1
    blitz::Array<int16_t,2> zicetop1m(JM1m, IM1m);
    {fortran::UnformattedInput fin(files.locate("ZETOPO1.NCEI"), Endian::BIG);

        // Continental cells north of 78N are entirely glacial ice.
        // (Except for Greenland, if that is to be done on a separate ice sheet)
        {
            blitz::Array<int16_t,2> focean1m(JM1m, IM1m);
            fortran::read(fin) >> titlei >> focean1m >> fortran::endr;

            // Continental cells north of 78N are entirely glacial ice.
            // (but ignore Greenland)
            for (int j1m=JM1m*14/15; j1m < JM1m; ++j1m) {
            for (int i1m=0; i1m < IM1m; ++i1m) {
                if (focean1m(j1m,i1m) == 0) {
                    // Convert indexing from 1minute --> 1degree grid
                    int const i1 = i1m / 60;
                    int const j1 = j1m / 60;

                    if (include_greenland || fcont1(j1,i1) != 2)
                        fgice1m(j1m, i1m) = 1;
                }
            }}

            {NcIO ncout(fname_out, 'w');
                ncio_blitz(ncout, focean1m, "FOCEAN1m", "int16", {});
            }

        }


        // Antarctica is supplied by ETOPO1
        {
            blitz::Array<int16_t,2> zsolid1m(JM1m, IM1m);

            fortran::read(fin) >> titlei >> zicetop1m >> fortran::endr;
            fortran::read(fin) >> titlei >> zsolid1m >> fortran::endr;

            {NcIO ncout(fname_out, 'a');
                ncio_blitz(ncout, zicetop1m, "ZICETOP1m", "int16", {});
                ncio_blitz(ncout, zsolid1m, "ZSOLG1m", "int16", {});
            }


            // Use ETOPO1 for Southern Hemisphere Ice
            for (int j1m=0; j1m < JM1m/2; ++j1m) {
            for (int i1m=0; i1m < IM1m; ++i1m) {
                if (zicetop1m(j1m,i1m) != zsolid1m(j1m,i1m))
                    fgice1m(j1m, i1m) = 1;
            }}

            // (Maybe) use EOTOPO1 for Greenland too
            if (include_greenland) {
                for (int j1m=JM1m/2; j1m < JM1m; ++j1m) {
                for (int i1m=0; i1m < IM1m; ++i1m) {
                    if (zicetopI1m(j1m,i1m) != zsolid1m(j1m,i1m))
                        fgice1m(j1m, i1m) = 1;
                }}
            }
        }
    }

    // Add northern-hemisphere non-Greenland ice specified in ZNGDC1
    // This must be downscaled to ETOPO1 grid
    std::vector<std::tuple<int16_t,int,int>> cells;
    for (int j1=JM1/2; j1 < JM1; ++j1) {
    for (int i1=0; i1 < IM1; ++i1) {
        if (is_greenland1(j1,i1)) continue;


        // Assemble cells in this gridcell, sorted by
        // descending elevation (increasing negative elevation)
        for (j1m=j1*60; j1m < (j1+1)*60; ++j1m) {
        for (i1m=i1*60; i1m < (i1+1)*60; ++i1m) {
            cells.push_back(std::make_tuple(
                -zicetopIs(j1m,i1m), j1m, i1m));
        }}
        std::sort(cells.begin(), cells.end());

        // Snow-covered area for this cell in g1x1
        double snow1 = fgice1(j1, i1);
        double remain1m = snow1 * g1x1_dxyp(j1);
        for (auto ii=cells.begin(); ii != cells.end(); ++cells) {
            int16_t const elev1m = -std::get<0>(*ii);
            int const j1m(std::get<1>(*ii));
            int const i1m(std::get<2>(*ii));

            double const area1m = g1mx1m_dxyp(j1m);
            if (area1m <= remain1m) {
                remain1m = -area1m;
                fgice1m(j1m, i1m) = 1;
            } else {
                // We're done; round the last grid cell on g1mx1m
                if (remain1m >= area1m * .5)
                    fgice1m(j1m, i1m) = 1;
                break;
            }
        }
    }}


    // Save fgice1m
    {NcIO ncout(fname_out, 'a');
        ncio_blitz(ncout, fgice1m, "FGICE1m" "int16", {});
    }
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
            "ofname", "Name of output file", true, "", "output filename", cmd);
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
    etopo1_ice(FileLocator("MODELE_FILE_PATH"), args.greenland);
}


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

