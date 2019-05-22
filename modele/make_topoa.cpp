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

        // Parse the argv array.
        cmd.parse( argc, argv );

        topoo_fname = topoo_a.getValue();
//        global_ec_fname = global_ec_a.getValue();
        global_ec_mm_fname = global_ec_mm_a.getValue();
//        elevmask_fname = elevmask_a.getValue();
//            elevI_vname = elevI_vname_a.getValue();
//            fgiceI_vname = fgiceI_vname_a.getValue();

        topoa_fname = topoa_a.getValue();
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

    // =========== Read grids, etc.
    global_ec::Metadata meta;
    {NcIO ncio(args.global_ec_mm_fname, 'r');
        meta.ncio(ncio);
    }

    HntrSpec hspecO;
    {NcIO ncio(args.topoo_fname, 'r');
        hspecO.ncio(ncio, "hspec");
    }

    // Create the AvO regridder (but not matrix)
    Hntr hntr_AvO(17.17, meta.hspecA, hspecO);

    // ============= Convert TOPOO to TOPOA
    std::map<std::string, std::unique_ptr<blitz::Array<double,2>>> varsA;

    ibmisc::ArrayBundle<double,2> topoo;
    auto &foceanOm(out.add("FOCEAN", {
        "description", "0 or 1, Bering Strait 1 cell wide",
        "units", "1",
        "source", "GISS 1Qx1",
    }));
    auto &flakeOm(out.add("FLAKE", {
        "description", "Lake Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &fgrndOm(out.add("FGRND", {
        "description", "Ground Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &fgiceOm(out.add("FGICE", {
        "description", "Glacial Ice Surface Fraction",
        "units", "0:1",
        "sources", "GISS 1Qx1",
    }));
    auto &zatmoOm(out.add("ZATMO", {
        "description", "Atmospheric Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &zlakeOm(out.add("ZLAKE", {
        "description", "Lake Surface Topography",
        "units", "m",
        "sources", "ETOPO2 1Qx1",
    }));
    auto &zicetopOm(out.add("ZICETOP", {
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
    auto &foceanOp(out.add("FOCEANF", {
        "description", "0 or 1, Bering Strait 1 cell wide",
        "units", "1",
        "source", "GISS 1Qx1",
    }));


    topoa.allocate(true);    // check=true
    auto &foceanA(topoa.array("focean"));
    auto &flakeA(topoa.array("flake"));
    auto &fgrndA(topoa.array("fgrnd"));
    auto &fgiceA(topoa.array("fgice"));
    auto &zatmoA(topoa.array("zatmo"));
    auto &zlakeA(topoa.array("zlake"));
    auto &zicetopA(topoa.array("zicetop"));


    // Read TOPOO input
    {NcIO topoo_nc(args.topoo_fname, 'r');

        // Read from topoO file, and allocate resulting arrays.
        topoo.ncio_alloc(topoo_nc, {}, "", "double",
            get_or_add_dims(topoo_nc, {"jm", "im"}));
    }

    // Create gcmA

    // Regrid to A grid, etc.
    make_topoA(&*gcmA,
        foceanOp, foceanOm, flakeOm, fgrndOm, fgiceOm, zatmoOm, zlakeOm, zicetopOm,
        EOpvAOp, paramsA, dimeEOp, dimAOp, eq_rad,
        foceanA, flakeA, fgrndA, fgiceA, zatmoA, zlakeA, zicetopA,
        fhc, elev, underice);

    // Write extended TOPO file
    {NcIO topoa_nc(args.topoa_fname, 'w');

        // Read from topoO file, and allocate resulting arrays.
        topoa.ncio(topoa_nc, {}, "", "double",
            get_or_add_dims(topoa_nc, {"jm", "im"}));
    }












        // Set up weight vectors for each regrid
        blitz::Array<double,2> fgiceO(hspecO.jm, hspecO.im);
        ncio_blitz(topoo_nc, fgiceO, "FGICE", "", {});
        blitz::Array<double, 2> WTO(const_array(blitz::shape(hspecO.jm,hspecO.im), 1.0));
        std::vector<blitz::Array<double,2> *> const weightsO
            {&WTO, &WTO, &WTO, &WTO, &WTO, &WTO, &fgiceO};

        // Define output variables
        auto dimsA(get_or_add_dims(topoa_nc, {"jm", "im"}, {meta.hspecA.jm, meta.hspecA.im}));

        // Create the AvO regridder (but not matrix)
        Hntr hntr_AvO(17.17, meta.hspecA, hspecO);

        // Regrid TOPOO variables and save to TOPOA
        for (int i=0; i<itopo_vars.size(); ++i) {
            printf("Reading %s (O)\n", itopo_vars[i].c_str());

            // Read on O grid
            blitz::Array<double,2> valO2(hspecO.jm, hspecO.im);
            NcVar ncv(ncio_blitz(topoo_nc, valO2, itopo_vars[i], "", {}));

            // Allocate variable to hold it
            std::unique_ptr<blitz::Array<double,2>> varAp(
                new blitz::Array<double,2>(meta.hspecA.jm, meta.hspecA.im));
            auto &varA(*varAp);
            varsA.insert(std::make_pair(otopo_vars[i], std::move(varAp)));

            // Regrid to A grid
            hntr_AvO.regrid(*weightsO[i], valO2, varA);
        }

        for (int i=0; i<itopo_vars.size(); ++i) {
            printf("Writing %s (A)\n", otopo_vars[i].c_str());

            std::vector<std::pair<std::string, NcAttValue>> atts;
            NcVar ncv_in(topoo_nc.nc->getVar(itopo_vars[i]));
            get_or_put_all_atts(ncv_in, topoo_nc.rw, atts);


            // Write it out, transferring attributes from ori ginal variable
            auto &varA(*varsA.at(otopo_vars[i]));
            NcVar ncv_out(ncio_blitz(topoa_nc, varA, otopo_vars[i], "double", dimsA));
            get_or_put_all_atts(ncv_out, topoa_nc.rw, atts);
        }
    }

    // ================= Create fhc, elevE and underice
    int const nhc_icebin = meta.indexingE[2].extent;
    int const nhc_gcm = 1 + 1 + nhc_icebin;
    auto dimsE(get_or_add_dims(topoa_nc,
        {"nhc", "jm", "im"},
        {nhc_gcm, meta.hspecA.jm, meta.hspecA.im}));
    blitz::TinyVector<int,3> shapeE(nhc_gcm, meta.indexingE[1].extent, meta.indexingE[0].extent);
    blitz::TinyVector<int,2> shapeE2(shapeE(0), shapeE(1)*shapeE(2));

printf("shapeE = (%d, %d, %d)\n", shapeE(0), shapeE(1), shapeE(2));
printf("shapeE2 = (%d, %d)\n", shapeE2(0), shapeE2(1));


    blitz::Array<double,3> fhc(shapeE);
    blitz::Array<double,3> elevE(shapeE);
    blitz::Array<uint16_t,3> underice(shapeE);

    // Initialize
    fhc = 0;
    elevE = NaN;
    underice = 0;    // Elevation grid cell unused

    auto all(blitz::Range::all());
    int ec_base = 0;


    // ------------ Segment 0: Legacy
    printf("Segment 0: Legacy\n");
    auto &fgiceA(*varsA.at("fgice"));
    auto &zatmoA(*varsA.at("zatmo"));
    auto &zicetopA(*varsA.at("zicetop"));

    for (int j=0; j<meta.hspecA.jm; ++j) {
    for (int i=0; i<meta.hspecA.im; ++i) {
        if (fgiceA(j,i) > 0) {
            fhc(ec_base, j,i) = 1. * 1e-30;
            elevE(ec_base, j,i) = zatmoA(j,i);
            underice(ec_base, j,i) = UI_NOTHING;
        }
    }}
    ec_base += 1;

    // ------------ Segment 1: land part of sealand
    printf("Segment 1: seaLAND\n");
    for (int j=0; j<meta.hspecA.jm; ++j) {
    for (int i=0; i<meta.hspecA.im; ++i) {
        if (fgiceA(j,i) > 0) {
            fhc(ec_base, j,i) = 1e-30;
            elevE(ec_base, j,i) = zicetopA(j,i);
            underice(ec_base, j,i) = UI_NOTHING;
        }
    }}
    ec_base += 1;

    // ------------ Segment 2: Elevation Classes
    printf("Segment 2: Elevation Classes\n");
    auto fhcE2(reshape<double,3,2>(fhc, shapeE2));
    auto elevE2(reshape<double,3,2>(elevE, shapeE2));
    auto undericeE2(reshape<uint16_t,3,2>(underice, shapeE2));
    {NcIO ncio(files.locate(args.global_ec_mm_fname), 'r');
        auto AvE(linear::nc_read_weighted(ncio.nc, "AvE"));

        // Uncompress AvE
        std::array<blitz::Array<int,1>,2> indices;
        blitz::Array<double,1> values;
        AvE->to_coo(indices[0], indices[1], values);

        // Scan AvE to get fhc
        int const nnz(AvE->nnz());
        std::array<int,2> iTuple;
            int &iA2(iTuple[0]);
            int &ihc(iTuple[1]);

        for (int ihc=0; ihc<meta.hcdefs.size(); ++ihc) {
            printf("elevE[ihc=%d] = %g\n", ec_base+ihc, meta.hcdefs[ihc]);
        }

        for (int i=0; i<nnz; ++i) {
            auto const iA(indices[0](i));
            auto const iE(indices[1](i));

            // iE must be contained within cell iA (local propety of matrix)
            meta.indexingHC.index_to_tuple(&iTuple[0], iE);

            if (iA2 != iA) (*icebin_error)(-1,
                "Matrix is non-local: iA=%ld, iE=%ld, iA2=%ld",
                (long)iA, (long)iE, (long)iA2);

            if (ihc < 0 || ihc >= nhc_icebin) (*icebin_error)(-1,
                "ihc out of range [0,%d): %d", nhc_icebin, ihc);

            if (iA < 0 || iA >= shapeE2(1)) (*icebin_error)(-1,
                "iA out of range [0,%d): %d", shapeE2(1), iA);

if (values(i) < 0 || values(i) >= 1.) printf("AvE(%d, ihc=%d) = %g\n", iA, ihc, values(i));


            fhcE2(ec_base+ihc,iA) += values(i);
            undericeE2(ec_base+ihc,iA) = UI_NOTHING;    // No IceBin coupling here
            elevE2(ec_base+ihc,iA) = meta.hcdefs[ihc];
        }
    }

    // --------------- Write it out
    NcVar info(get_or_add_var(topoa_nc, "info", "int", {}));
    std::string sval = "legacy,land,ec";
    get_or_put_att(info, topoa_nc.rw, "segments", sval);

    NcVar ncv;
    ncio_blitz(topoa_nc, fhc, "fhc", "double", dimsE);
    ncio_blitz(topoa_nc, elevE, "elevE", "double", dimsE);
    ncio_blitz(topoa_nc, underice, "underice", "double", dimsE);

    topoa_nc.flush();
    topoa_nc.close();

    return 0;
}
