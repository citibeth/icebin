#include <algorithm>
#include <string>
#include <iostream>
#include <tclap/CmdLine.h>
#include <ibmisc/memory.hpp>
#include <ibmisc/blitz.hpp>
#include <ibmisc/filesystem.hpp>
#include <ibmisc/linear/linear.hpp>
#include <everytrace.h>

#include <icebin/modele/grids.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/global_ec.hpp>

using namespace ibmisc;
using namespace icebin;
using namespace icebin::modele;

static double const NaN = std::numeric_limits<double>::quiet_NaN();

struct ParseArgs {
    std::string topoo_fname;
//    std::string global_ec_fname;
    std::string global_ec_mm_fname;
// These variables are vestigal
//    std::string elevmask_fname;
//        std::string elevI_vname, fgiceI_vname;

    std::string topoa_fname;
    double eq_rad;

    ParseArgs(int argc, char **argv);
};

ParseArgs::ParseArgs(int argc, char **argv)
{
    // Wrap everything in a try block.  Do this every time, 
    // because exceptions will be thrown for problems.
    try {  
        // Define the command line object, and insert a message
        // that describes the program. The "Command description message" 
        // is printed last in the help text. The second argument is the 
        // delimiter (usually space) and the last one is the version number. 
        // The CmdLine object parses the argv array based on the Arg objects
        // that it contains. 
        TCLAP::CmdLine cmd("Command description message", ' ', "<no-version>");

        TCLAP::ValueArg<std::string> topoo_a("a", "topoo",
            "TOPOO file, writen by make_topoo",
            false, "topoo.nc", "topoo file", cmd);

#if 0
        TCLAP::ValueArg<std::string> global_ec_a("b", "global_ec",
            "Elevation Class Matrix file (NOT mismatched)",
            false, "global_ec.nc", "topoo file", cmd);
#endif

        TCLAP::ValueArg<std::string> global_ec_mm_a("c", "global_ec_mm",
            "Elevation Class Matrix file (mismatched)",
            false, "global_ec_mm.nc", "topoo file", cmd);

        TCLAP::ValueArg<std::string> elevmask_a("d", "elevmask",
            "Source file for FGICE1m and ZICETOP1m",
            false, "FGICE1m", "focean var name", cmd);

            TCLAP::ValueArg<std::string> elevI_vname_a("e", "elev",
                "Name of NetCDF variable containing elevation",
                false, "ZICETOP1m", "focean var name", cmd);

            TCLAP::ValueArg<std::string> fgiceI_vname_a("f", "mask",
                "Name of NetCDF variable containing ice mask (1 where there is ice)",
                false, "FGICE1m", "mask var name", cmd);

        TCLAP::ValueArg<std::string> topoa_a("o", "topoa",
            "OUT: TOPO file on Atmosphere grid",
            true, "global_ec.nc", "topoa file", cmd);

        TCLAP::ValueArg<double> eq_rad_a("R", "radius",
            "Radius of the earth",
            false, modele::EQ_RAD, "earth radius", cmd);


        // Parse the argv array.
        cmd.parse( argc, argv );

        topoo_fname = topoo_a.getValue();
//        global_ec_fname = global_ec_a.getValue();
        global_ec_mm_fname = global_ec_mm_a.getValue();
//        elevmask_fname = elevmask_a.getValue();
//            elevI_vname = elevI_vname_a.getValue();
//            fgiceI_vname = fgiceI_vname_a.getValue();

        topoa_fname = topoa_a.getValue();
        eq_rad = eq_rad_a.getValue();
    } catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        exit(1);
    }
}




int main(int argc, char **argv)
{
    everytrace_init();
    ParseArgs args(argc, argv);
    EnvSearchPath files("MODELE_FILE_PATH");

    // Open output file
    NcIO topoa_nc(args.topoa_fname, 'w');


    // =========== Read metadata and EOpvAOp matrix
    global_ec::Metadata metaO;
    std::unique_ptr<EigenSparseMatrixT> EOpvAOp;
    {NcIO ncio(args.global_ecO, 'r');
    linear::Weighted_Compressed EOpvAOp_s;    // sparse indexing
        metaO.ncio(ncio);

        EOpvAOp_s.ncio(ncio, "EvA");
        EOpvAOp.reset(new EigenSparseMatrixT(to_eigen_M(EOpvAOp_s)));
    }

    HntrSpec &hspecO(metaO.hspecA);
    HntrSpec hspecA(make_hntrA(hspecO));
    Indexing &indexingHCO(metaO.indexingHC);
    Indexing indexingHCA({"A", "HC"}, {0,0}, {hspecA.size(), indexingHCO[1].extent}, {1,0}),

    // Create the AvO regridder (but not matrix)
    Hntr hntr_AvO(17.17, meta.hspecA, hspecO);

    // ============= Define input variables
    ibmisc::ArrayBundle<double,2> topoo;
    auto &foceanOm(topoo.add("FOCEAN", {
        "description", "0 or 1, Bering Strait 1 cell wide",
        "units", "1",
        "source", "GISS 1Qx1",
    }));
    auto &flakeOm(topoo.add("FLAKE", {
        "description", "Lake Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &fgrndOm(topoo.add("FGRND", {
        "description", "Ground Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &fgiceOm(topoo.add("FGICE", {
        "description", "Glacial Ice Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &zatmoOm(topoo.add("ZATMO", {
        "description", "Atmospheric Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &zlakeOm(topoo.add("ZLAKE", {
        "description", "Lake Surface Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &zicetopOm(topoo.add("ZICETOP", {
        "description", "Atmospheric Topography (Ice-Covered Regions Only)",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));


    // Copy TOPOO to TOPOA spec, downcasing variable names
    ibmisc::ArrayBundle<double,2> topoa(topoo);
    for (auto &d : topoa.data) {
        std::string &name(d.meta.name);
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);

        // Set shape
        d.meta.shape = std::array<int,2>{meta.hspecA.jm, mta.hspecA.im};
    }

    // Add FOCEANF, which is in TOPOO but not TOPOA.
    auto &foceanOp(topoo.add("FOCEANF", {
        "description", "0 or 1, Bering Strait 1 cell wide",
        "units", "1",
        "source", "GISS 1Qx1",
    }));


    // Read TOPOO input
    {NcIO topoo_nc(args.topoo_fname, 'r');

        // Read from topoO file, and allocate resulting arrays.
        topoo.ncio_alloc(topoo_nc, {}, "", "double",
            get_or_add_dims(topoo_nc, {"jm", "im"}));
    }


    // ------------- Allocate topoa arrays for output
    topoa.allocate(true);    // check=true
    auto &foceanA(topoa.array("focean"));
    auto &flakeA(topoa.array("flake"));
    auto &fgrndA(topoa.array("fgrnd"));
    auto &fgiceA(topoa.array("fgice"));
    auto &zatmoA(topoa.array("zatmo"));
    auto &zlakeA(topoa.array("zlake"));
    auto &zicetopA(topoa.array("zicetop"));

    // --------------- Allocate 3D arrays to go in TOPOA file
    ibmisc::ArrayBundle<double,3> topoa3;
    auto &fhc(topoa3.add("fhc", {
        "description", "fraction of ice-covered area for each elevation class",
        "units", "1"
    }));
    auto &elevE(topoa3.add("elevE", {
        "description", "Elevation of each elevation class",
        "units", "1"
    }));
    ibmisc::ArrayBundle<double,3> topoa3_i;
    auto &underice(topoa3_i.add("underice", {
        "description", "Model below the show/firn (UI_UNUSED=0, UI_ICEBIN=1, UI_NOTHING=2)"
    }));


    int const nhc_gcm = meta.hcdefs.
    std::array<int,3> shape3(nhc_gcm, metaO.indexingA[1].extent, meta.indexingA[0].extent);
    topoa3.allocate(shape3);
    topoa3_i.allocate(shape3);


    // ---------------- Create TOPOA in memory
    std::vector<uint16_t> underice;
    for (size_t i=0; i<meta.hcdefs.size(); ++i) underice.push_back(UI_NOTHING);
    make_topoA(
        foceanOp, foceanOm, flakeOm, fgrndOm, fgiceOm, zatmoOm, zlakeOm, zicetopOm,
        hspecO, hspecA, indexingHCO, indexingHCA, meta.hcdefs, underice,
        eq_rad, EOpvAOp, dimeEOp, dimAOp,
        foceanA, flakeA, fgrndA, fgiceA, zatmoA, zlakeA, zicetopA,
        fhc, elev, underice);

    // Write extended TOPOA file
    {NcIO topoa_nc(args.topoa_fname, 'w');

        NcVar info(get_or_add_var(topoa_nc, "info", "int", {}));
        std::string sval = "ec,land";
        get_or_put_att(info, topoa_nc.rw, "segments", sval);

        // Write topoA arrays
        auto jm_im(get_or_add_dims(topoa_nc, {"jm", "im"}));
        topoa.ncio(topoa_nc, {}, "", "double", jm_im);
        topoa3.ncio(topoa_nc, {}, "", "double", jm_im);
        topoa3_i.ncio(topoa_nc, {}, "", "ushort", jm_im);
    }

    return 0;
}
