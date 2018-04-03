#include <icebin/modele/grids.hpp>
#include <icebin/modele/hntr.hpp>


TOPOO
GLOBAL_EC
GLOBAL_EC_MISMATCHED
etopo1_ice_g1m.nc
OUT: TOPOA

struct ParseArgs {
    std::string topoo_fname;
    std::string global_ec_fname;
    std::string global_ec_mm_fname;
    std::string elevmakI_fname;
        std::string elevI_vname, fgiceI_vname;

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

        TCLAP::ValueArg<std::string> ecmatrix_a("b", "ecmatrix",
            "Elevation Class Matrix file (NOT mismatched)",
            false, "global_ec.nc", "topoo file", cmd);

        TCLAP::ValueArg<std::string> ecmatrix_mm_a("c", "ecmatrix_mm",
            "Elevation Class Matrix file (mismatched)",
            false, "global_ec_mm.nc", "topoo file", cmd);

        TCLAP::ValueArg<std::string> elevmask_a("d", "elevmask",
            "Source for FGICE1m and ZICETOP1m",
            false, "FGICE1m", "focean var name", cmd);

            TCLAP::ValueArg<std::string> elevI_vname_a("e", "elev",
                "Name of NetCDF variable containing elevation"
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
        ecmatrix_fname = ecmatrix_a.getValue();
        ecmatrix_mm_fname = ecmatrix_mm_a.getValue();
        elevmaskI_fname = elevmask_a.getValue();
            elevI_vname = elevI_vname_a.getValue();
            fgiceI_vname = fgiceI_vname_a.getValue();

        topoa_fname = topoa_a.getValue();
    } catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        exit(1);
    }
}


static std::vector<std::string> const itopo_vars
    {"FOCEAN", "FLAKE", "FGRND", "FGICE", "ZATMO", "ZLAKE", "ZICETOP"};
static std::vector<std::string> const otopo_vars
    {"focean", "flake", "fgrnd", "fgice", "zatmo", "hlake", "zicetop"};



void make_topoa(ParseArgs args)
{

    // Open output file
    NcIO topoa_nc(args.topoa_fname, 'w');

    // ============= Convert TOPOO to TOPOA
    blitz::Array<double,2> foceanA, flakeA, fgrndA, fgiceA, zatmoA, hlakeA;
    std::vector<blitz::Array<double,2> *> const varsA
        {&foceanA, &flakeA, &fgrndA, &fgiceA, &zatmoA, &hlakeA, &zicetopA};
    {NcIO topoo_nc(TOPOO_nc, 'r');

        // Set up weight vectors for each regrid
        blitz::Array<double,2> fgiceO(hspecO.jm, hspecO.im);
        ncio_blitz(topoo_nc, fgiceO, "FGICE", "");
        blitz::Array<double, 2> WTO(const_array(hspecO.jm,hspecO.im, 1.0));
        std::vector<blitz::Array<double,2> *> const weightsO
            {&WTO, &WTO, &WTO, &WTO, &WTO, &WTO, &fgiceO};

        // Define output variables
        auto adims(get_or_add_dims(topoa_nc, {"jm", "im"}, {hspecA.jm, hspecA.im}));

        // Create the AvO matrix (Eigen format)
        Hntr hntr_AvO(17.17, args.hspecA, args.hspecO);
#if 0
        TupleListT<2> AvO_tp;
        hntrAvO.scaled_regrid_matrix(spsparse::accum::ref(AvO_tp));
        EigenSparseMatrixT AvO_e(hntrA.size(), hntrO.size());
        AvO_e.setFromTriplets(AvO_tp.begin(), AvO_tp.end());
#endif

        // Regrid TOPOO variables and save to TOPOA
        for (int i=0; i<itopo_vars.size(); ++i) {
            // Read on O grid
            blitz::Array<double,2> valO2(hspecO.jm, hspecO.im);
            std::vector<std::pair<std::string, NcAttValue>> atts;
            get_or_put_all_atts(
                ncio_blitz(topoo_nc, valO2, itopo_vars[i], ""),
                topoo_nc.rw, atts);
            blitz::Array<double,1> valO(reshape1(valO2));

            // Regrid to A grid
            auto &valA2(topoa_nc.tmp.make<blitz::Array<double,2>>(hspecA.jm, hspecA.im));
            blitz::Array<double,1> valA(reshape1(valA2));
            if (otopo_vars[i] == "zicetop") {
                hntrAvO.regrid(*weightsO[i], valO, valA);
            }
//                map_eigen_colvector(valA) = AvO_e * map_eigen_colvector(valO);

            // Save for other computations
            varsA[i]->reference(valA2);
        }

        // ZICETOP only exists in ice-covered regions, so it must be weighted differently.
        // (This is a cheap trick; we really should use hntr.regrid(), with
        // WT=FGICE.  But as long as ZICETOP=0 over non-ice, this trick will work
        if (otopo_vars[i] == "zicetop") {
            for (int jj=0; jj<zicetopA.extent(0); ++jj) {
            for (int ii=0; ii<zicetopA.extent(1); ++ii) {
                if (fgiceA(jj,ii) > 1.e-8) {
                    zicetopA(jj,ii) = zicetopA(jj,ii) / fgiceA(jj,ii);
                }
            }}
        }

        for (int i=0; i<itopo_vars.size(); ++i) {
            // Write it out, transferring attributes from ori ginal variable
            get_or_put_all_atts(
                ncio_blitz(topoa_nc, *varsA[i], otopo_vars[i], "double", adims),
                topoa_nc.rw, attss);
        }

    }

    // ================= Create fhc, elevE and underice
    int const nhc_icebin = indexingE[2].extent;
    int const nhc_gcm = 1 + 1 + nhc_icebin;
    blitz::TinyVector<int,3> shapeE(indexingE[2].extent, indexingE[1].extent, indexingE[0].extent);
    blitz::TinyVector<int,3> shapeE2(shapeE(0), shapeE(1)*shapeE(2));


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
    for (int j=0; j<hspecA.jm; ++j) {
    for (int i=0; i<hspecA.im; ++i) {
        if (fgiceA(j,i) > 0) {
            fhc(ec_base, j,i) = 1. * 1e-30;
            elevE(ec_base, j,i) = zatmoA(j,i);
            underice(ec_base, j,i) = UI_NOTHING;
        }
        ec_base += 1;
    }}

    // ------------ Segment 1: land part of sealand
    for (int j=0; j<hspecA.jm; ++j) {
    for (int i=0; i<hspecA.im; ++i) {
        if (fgiceA(j,i) > 0) {
            fhc(ec_base, j,i) = 1e-30;
            elevE(ec_base, j,i) = zicetopA(j,i);
            underice(ec_base, j,i) = UI_NOTHING;
        }
        ec_base += 1;
    }}


    // ------------ Segment 2: Elevation Classes
    auto fhcE2(reshape<double,3,2>(fhc, shapeE2);
    {NcIO ncio(files.locate(args.global_ec_mm), 'r');
        auto AvE(nc_read_weighted(ncio.nc, "AvE"));

        // Uncompress AvE
        std::array<blitz::Array<int,1>,2> indices;
        blitz::Array<double,2> values;
        AvE->to_coo(indices[0], indices[1], values);

        // Scan AvE to get fhc
        int const nnz(AvE->nnz());
        std::array<int,2> iTuple;
            int &iA2(iTuple[0]);
            int &ihc(iTuple[1]);
        for (int i=0; i<nnz; ++i) {
            auto const iA(indices[0](i));
            auto const iE(indices[0](i));

            // iE must be contained within cell iA (local propety of matrix)
            indexingHC.index_to_tuple(&iTuple[0], iE);

            if (iA2 != iA) (*icebin_error)(-1,
                "Matrix is non-local: iA=%ld, iE=%ld, iA2=%ld",
                (long)iA, (long)iE, (long)iA2);

            fhcE2(ec_base+ihc,iA) = values(i);
            underice(ec_base+ihc,iA) = UI_ICEBIN;
            elevE2(ec_base+ihc,iA) = meta.hcdefs[ihc];
            underice(ec_base, all, all) = UI_NOTHING;
        }
    }

    // --------------- Write it out
    {NcIO ncio(args.topoa_fname, 'w');
        NcVar info(get_or_add_var(ncio, "info", "int", {}));
        get_or_put_att(info, ncio.rw, "segments", "legacy,sealand,ec");



    }
}


