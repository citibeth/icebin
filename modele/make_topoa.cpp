#include <algorithm>
#include <string>
#include <iostream>
#include <tclap/CmdLine.h>
#include <ibmisc/memory.hpp>
#include <ibmisc/blitz.hpp>
#include <ibmisc/filesystem.hpp>
#include <ibmisc/linear/compressed.hpp>
#include <ibmisc/linear/eigen.hpp>
#include <everytrace.h>
#include <boost/filesystem.hpp>

#include <icebin/modele/grids.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/global_ec.hpp>
#include <icebin/modele/topo.hpp>

using namespace netCDF;
using namespace ibmisc;
using namespace icebin;
using namespace icebin::modele;

static double const NaN = std::numeric_limits<double>::quiet_NaN();

struct ParseArgs {
    std::string topoo_fname;
//    std::string global_ec_fname;
    std::string global_ecO_fname;
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

        TCLAP::ValueArg<std::string> global_ecO_a("c", "global_ecO",
            "Elevation Class Matrix file (elevation grid)",
            false, "global_ecO.nc", "matrix file", cmd);

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
        global_ecO_fname = global_ecO_a.getValue();
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

    // =========== Read metadata and EOpvAOp matrix
    std::unique_ptr<EigenSparseMatrixT> EOpvAOp;
    SparseSetT dimEOp, dimAOp;
    icebin::HntrSpec hspecO;
    ibmisc::Indexing indexingHC;
    std::vector<double> hcdefs;
    std::vector<int16_t> underice_hc;
    {NcIO ncio(args.global_ecO_fname, 'r');
    ZArray<int,double,2> EOpvAOp_s;

        hspecO.ncio(ncio, "hspecO");    // Actually ocean grid
        indexingHC.ncio(ncio, "indexingHC");
        ncio_vector(ncio, hcdefs, true, "hcdefs", "double", {});
        ncio_vector(ncio, underice_hc, true, "underice_hc", "short", {});  // Must be short for NetCDF3

        EOpvAOp_s.ncio(ncio, "EvO.M");
        EOpvAOp.reset(new EigenSparseMatrixT(
            to_eigen_M(EOpvAOp_s, {&dimEOp, &dimAOp})));
    }
    int const nhc_gcm = get_nhc_gcm(hcdefs.size());

    HntrSpec hspecA(make_hntrA(hspecO));
    Indexing &indexingHCO(indexingHC);
    Indexing indexingHCA({"A", "HC"}, {0,0}, {hspecA.size(), indexingHCO[1].extent}, {1,0});

    // Create the AvO regridder (but not matrix)
    Hntr hntr_AvO(17.17, hspecA, hspecO);


    // Read input file, and allocate output arrays, ready to save to output file.
    auto topoo(topoo_bundle(BundleOType::MAKEA, args.topoo_fname));
    TopoABundles topoa(topoo, hspecA, nhc_gcm);
    blitz::Array<int16_t,2> mergemaskOm(hspecO.jm, hspecO.im);
    {NcIO ncio(args.topoo_fname, 'r');
        ncio_blitz(ncio, mergemaskOm, "MERGEMASK", "short", {});
    }


    // --------- Fetch just-allocated arrays
    auto &foceanOm(topoo.array("FOCEAN"));
    auto &flakeOm(topoo.array("FLAKE"));
    auto &fgrndOm(topoo.array("FGRND"));
    auto &fgiceOm(topoo.array("FGICE"));
    auto &zatmoOm(topoo.array("ZATMO"));
    auto &zlakeOm(topoo.array("ZLAKE"));
    auto &zicetopOm(topoo.array("ZICETOP"));
    auto &foceanOp(topoo.array("FOCEANF"));

    auto &foceanA(topoa.a.array("focean"));
    auto &flakeA(topoa.a.array("flake"));
    auto &fgrndA(topoa.a.array("fgrnd"));
    auto &fgiceA(topoa.a.array("fgice"));
    auto &zatmoA(topoa.a.array("zatmo"));
    auto &hlakeA(topoa.a.array("hlake"));
    auto &zicetopA(topoa.a.array("zicetop"));
    auto &mergemaskA(topoa.a_i.array("mergemask"));

    auto &fhc(topoa.a3.array("fhc"));
    auto &elevE(topoa.a3.array("elevE"));

    auto &underice(topoa.a3_i.array("underice"));

    // ------------------- Create AAmvEAm, needed for TOPOA
    std::unique_ptr<ConstUniverseT> const_dimAOp(
        new ConstUniverseT({"dimEOp", "dimAOp"}, {&dimEOp, &dimAOp}));

    // Compute AAmvEAm --> fhc
    auto wAOp(sum(*EOpvAOp, 1, '+'));
    linear::Weighted_Tuple AAmvEAm(_compute_AAmvEAm(
        true, args.eq_rad,    // scale=true
        hspecO, hspecA, indexingHCO, indexingHCA,
        reshape1(foceanOp), reshape1(foceanOm),   // These don't change over course of a run
        *EOpvAOp, dimEOp, dimAOp, wAOp));

    const_dimAOp.reset();    // Check for any const violations


    // ---------------- Create TOPOA in memory
    std::vector<std::string> errors(make_topoA(
        foceanOm, flakeOm, fgrndOm, fgiceOm, zatmoOm, zlakeOm, zicetopOm, mergemaskOm,
        hspecO, hspecA, indexingHCA, hcdefs, underice_hc,
        AAmvEAm,
        foceanA, flakeA, fgrndA, fgiceA, zatmoA, hlakeA, zicetopA, mergemaskA,
        fhc, elevE, underice));

    // Print sanity check errors to STDERR
    for (std::string const &err : errors) fprintf(stderr, "ERROR: %s\n", err.c_str());

    // Write extended TOPOA file
    std::vector<double> lonc(hspecA.lonc());
    std::vector<double> latc(hspecA.latc());
    {NcIO topoa_nc(args.topoa_fname, 'w');
        NcVar info(get_or_add_var(topoa_nc, "info", "int", {}));
        std::string sval = "ec,land";
        get_or_put_att(info, topoa_nc.rw, "segments", sval);

        ncio_vector(topoa_nc, lonc, false, "lon", "double",
            get_or_add_dims(topoa_nc, {"im"}, {lonc.size()}));
        ncio_vector(topoa_nc, latc, false, "lat", "double",
            get_or_add_dims(topoa_nc, {"jm"}, {latc.size()}));


        // Write topoA arrays
        auto nhc_jm_im(get_or_add_dims(topoa_nc, {"nhc", "jm", "im"}, {nhc_gcm, hspecA.jm, hspecA.im}));
        topoa.ncio(topoa_nc, nhc_jm_im);
    }

    if (errors.size() > 0) return -1;
}
