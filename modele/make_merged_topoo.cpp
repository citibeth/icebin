#include <algorithm>
#include <string>
#include <iostream>

#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>

#include <everytrace.h>

#include <ibmisc/memory.hpp>
#include <ibmisc/blitz.hpp>
#include <ibmisc/filesystem.hpp>
#include <ibmisc/linear/compressed.hpp>
#include <spsparse/eigen.hpp>

#include <icebin/GCMRegridder.hpp>
#include <icebin/modele/grids.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/global_ec.hpp>
#include <icebin/modele/topo.hpp>
#include <icebin/modele/merge_topo.hpp>

using namespace netCDF;
using namespace ibmisc;
using namespace spsparse;
using namespace icebin;
using namespace icebin::modele;

static double const NaN = std::numeric_limits<double>::quiet_NaN();

/** Command-line program reads:

    a) TOPOO (TOPO on oncean grid) and EOpvAOp matrix generated for
       global ice, but with ice sheet ("local ice") removed.

    b) GCMRegridder data structure, on the Ocean grid, capable of
       providing the missing ice sheet, directly from a hi-res form.
       For example, obtained from a PISM stat file.

Produces: TOPOO and EOpvAOp in which the local ice has been merged
    into the global ice.  This will later be processed by
    make_topoa.cpp to produce ModelE input files on the Atmosphere
    grid.
*/
struct ParseArgs {
    /** Name of the TOPOO file generated from the base ice.  It should
    be MISSING the ice sheets that will be provided by gcmO_fname.
    NOTE: `_ng` means "no Greenland" i.e. one or more ice sheets has
          been removed. */
    std::string topoo_ng_fname;

    /** Name of file out which the base EvA matrix (for global ice)
    will be loaded.  It should be MISSING the ice sheets that will be
    provided by gcmO_fname.
    NOTE: `_ng` means "no Greenland" i.e. one or more ice sheets has
          been removed. */
    std::string global_ecO_ng_fname;

    /** Name of file out of which the GCMRegirdder (ocean grid) for
    the local ice sheets will be loaded. */
    std::string gcmO_fname;

    /** Name of the files out of which the elevmaskI for each ice
    sheet will be loaded.  They must be in the same order as found in
    gcmO_fname.  Each filename is form of <format>:<fname>, allowing
    this program to know how to load and interpret the information in
    the file.  Currently, the only format is `pism:`; in which
    state files written by PISM are read. */
    std::vector<std::string> elevmask_xfnames;

    /** Should elevation classes between global and local ice be merged?
    This is desired when running without two-way coupling. */
    bool squash_ec;

    /** Output filename; the merged TOPOO and merged EOpvAOp matrices
    are written to this file. */
    std::string topoo_merged_fname;

    /** Radius of the earth to use when needed. */
    double eq_rad;

    ParseArgs(int argc, char **argv);
};

ParseArgs::ParseArgs(int argc, char **argv)
{
    // Wrap everything in a try block.  Do this every time, 
    // because exceptions will be thrown for problems.
    try {  
        TCLAP::CmdLine cmd("Command description message", ' ', "<no-version>");


        // ---------------- Input Filesnames
        TCLAP::ValueArg<std::string> topoo_ng_a("i", "topoo",
            "Knockout (eg Greenland-free) TOPOO file, writen by make_topoo",
            false, "topoo_ng.nc", "knockout topoo file", cmd);

        TCLAP::ValueArg<std::string> global_ecO_ng_a("c", "global_ecO",
            "Knockout (eg Greenland-free) Elevation Class Matrix file (ocean grid)",
            false, "global_ecO_ng.nc", "knockout matrix file", cmd);

        TCLAP::ValueArg<std::string> gcmO_ng_a("g", "gcmO",
            "File containing the GCMRegridder representing all ice sheets to be merged in (ocean grid)",
            false, "gcmO.nc", "GCMRegridder description", cmd);

        TCLAP::ValueArg<bool> squash_ec_a("s", "squash_ec",
            "Merge elevation classes between global and local ice?",
            false, true, "bool", cmd);

        TCLAP::MultiArg<std::string> elevmask_a("e", "elevmask",
            "<Source file for ice sheet elevation and maks>",
            false, "[type]:[fname]", cmd);

        TCLAP::ValueArg<double> eq_rad_a("R", "radius",
            "Radius of the earth",
            false, modele::EQ_RAD, "earth radius", cmd);


        // --------------- Output filenames
        TCLAP::ValueArg<std::string> topoo_merged_a("o", "topoo_merged",
            "Merged TOPOO file to write",
            false, "topoo_merged.nc", "output topoo file", cmd);

        // Parse the argv array.
        cmd.parse( argc, argv );

        // Extract values from TCLAP data structures.
        topoo_ng_fname = topoo_ng_a.getValue();
        global_ecO_ng_fname = global_ecO_ng_a.getValue();
        gcmO_fname = gcmO_ng_a.getValue();
        squash_ec = squash_ec_a.getValue();
        elevmask_xfnames = elevmask_a.getValue();
        topoo_merged_fname = topoo_merged_a.getValue();
    } catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        exit(1);
    }
}


int main(int argc, char **argv)
{
    everytrace_init();
    ParseArgs args(argc, argv);

    // ============= Define input/output  variables


    // ================================== Read Input Files

    // Read metadata and global EOpvAOp matrix (from output of global_ec.cpp)
    global_ec::Metadata metaO;
    ibmisc::ZArray<int,double,2> EOpvAOp_ng;
    {NcIO ncio(args.global_ecO_ng_fname, 'r');
        metaO.ncio(ncio);
        EOpvAOp_ng.ncio(ncio, "EvO.M");
    }
    HntrSpec &hspecO(metaO.hspecA);
    // HntrSpec hspecA(make_hntrA(hspecO));
    // Indexing &indexingHCO(metaO.indexingHC);
    // Indexing indexingHCA({"A", "HC"}, {0,0}, {hspecA.size(), indexingHCO[1].extent}, {1,0});

    // Read TOPOO input (global ice)
    ibmisc::ArrayBundle<double,2> topoo(topoo_bundle(BundleOType::MERGEO));
    blitz::Array<double,2> zsgloO;
        std::string zsgloO_units, zsgloO_sources;
    {NcIO topoo_nc(args.topoo_ng_fname, 'r');

        // Read from topoO file, and allocate resulting arrays.
        topoo.ncio_alloc(topoo_nc, {}, "", "double",
            get_or_add_dims(topoo_nc, {"jm", "im"}, {hspecO.jm, hspecO.im}));

        auto ncvar(ncio_blitz_alloc(topoo_nc, zsgloO, "ZSGLO", "double"));
            get_att(ncvar, "units").getValues(zsgloO_units);
            get_att(ncvar, "sources").getValues(zsgloO_sources);
    }
    auto &foceanOp(topoo.array("FOCEANF"));
    auto &fgiceOp(topoo.array("FGICEF"));
    auto &zatmoOp(topoo.array("ZATMOF"));
    auto &foceanOm(topoo.array("FOCEAN"));
    auto &flakeOm(topoo.array("FLAKE"));
    auto &fgrndOm(topoo.array("FGRND"));
    auto &fgiceOm(topoo.array("FGICE"));
    auto &zatmoOm(topoo.array("ZATMO"));
    auto &zlakeOm(topoo.array("ZLAKE"));
    auto &zicetopO(topoo.array("ZICETOP"));
    auto &zland_minO(topoo.array("ZLAND_MIN"));
    auto &zland_maxO(topoo.array("ZLAND_MAX"));

    // This is created in the merge.
    blitz::Array<int16_t,2> mergemaskO(hspecO.jm,hspecO.im);

    // Read the GCMRegridder
    GCMRegridder_Standard gcmO;
    {NcIO gcmO_nc(args.gcmO_fname, 'r');
        gcmO.ncio(gcmO_nc, "m");
    }

    // Read per-ice sheet elevmasks (for land+ice and ice only)
    std::vector<blitz::Array<double,1>> emI_lands, emI_ices;
    for (auto const &xfname : args.elevmask_xfnames) {
        blitz::Array<double,1> emI_land, emI_ice;
        read_elevmask(xfname, emI_land, emI_ice);

        // Store results
        emI_lands.push_back(emI_land);
        emI_ices.push_back(emI_ice);
    }

    std::vector<std::string> errors;

    // We need correctA=true here to get FOCEANF, etc.
    merge_topoO(
        foceanOp, fgiceOp, zatmoOp,
        foceanOm, flakeOm, fgrndOm, fgiceOm, zatmoOm, zicetopO,
        zland_minO, zland_maxO,
        mergemaskO, &gcmO,
        RegridParams(false, true, {0.,0.,0.}),  // (scale, correctA, sigma)
        emI_lands, emI_ices, args.eq_rad, errors);

    // Compute ZOCEAN by truncating 
    // This never changes in a model run (since the ocean is fixed)
    // Therefore, it doesn't need to be part of merge_topoO()
    blitz::Array<double,2> zoceanOm(hspecO.jm, hspecO.im);
    for (int j=0; j<hspecO.jm; ++j) {
    for (int i=0; i<hspecO.im; ++i) {
        zoceanOm(j,i) = foceanOm(j,i) == 1. ? std::max(0.,-zsgloO(j,i)) : 0.;
    }}



    SparseSetT dimAOp;
    EOpvAOpResult eam(compute_EOpvAOp_merged(
        dimAOp, EOpvAOp_ng,
        RegridParams(false, false, {0.,0.,0.}),  // (scale, correctA, sigma)
        &gcmO, args.eq_rad, emI_ices,
        true, true,    // use_global_ice=t, use_local_ice=t
        metaO.hcdefs, metaO.indexingHC, args.squash_ec, errors));


    // Print sanity check errors to STDERR
    for (std::string const &err : errors) fprintf(stderr, "ERROR: %s\n", err.c_str());

    // ================== Write output
    // Write all inputs to a single output file
    ZArray<int,double,2> EOpvAOp_c({eam.dimEOp.sparse_extent(), dimAOp.sparse_extent()});
    std::vector<double> lonc(metaO.hspecA.lonc());
    std::vector<double> latc(metaO.hspecA.latc());
    {NcIO ncio(args.topoo_merged_fname, 'w');

        // Write Ocean grid metadata
        metaO.hspecA.ncio(ncio, "hspecO");    // Actually ocean grid

        ncio_vector(ncio, lonc, false, "lon", "double",
            get_or_add_dims(ncio, {"im"}, {lonc.size()}));
        ncio_vector(ncio, latc, false, "lat", "double",
            get_or_add_dims(ncio, {"jm"}, {latc.size()}));

        eam.indexingHC.ncio(ncio, "indexingHC");
        auto xxdims(get_or_add_dims(ncio, {"nhc"}, {eam.hcdefs.size()}));
        ncio_vector(ncio, eam.hcdefs, true, "hcdefs", "double", xxdims);
        ncio_vector(ncio, eam.underice_hc, true, "underice_hc", "short", xxdims);  // Must be short for NetCDF3

        // Compress and Write EOpvAOp; our merged EOpvAOp needs to be
        // in the same (compressed) format as the original base
        // EOpvAOp that we read.
        {auto EOpvAOp_a(EOpvAOp_c.accum());
            for (auto ii=begin(*eam.EOpvAOp); ii != end(*eam.EOpvAOp); ++ii) {
                EOpvAOp_a.add({
                    eam.dimEOp.to_sparse(ii->index(0)),
                    dimAOp.to_sparse(ii->index(1))},
                    ii->value());
            }
        }    // Flush compression on ~EOpvAOp_a()

        // We write just the main matrix; but not the other things involved in
        // linear::Weighted_Compressed.
        EOpvAOp_c.ncio(ncio, "EvO.M");

        // Write out all the TOPOO items
        auto dims(get_or_add_dims(ncio, {"jm", "im"}, {hspecO.jm, hspecO.im}));
        topoo.ncio(ncio, {}, "", "double", dims);
        auto ncvar(ncio_blitz(ncio, zoceanOm, "ZOCEAN", "double", dims));
            ncvar.putAtt("description", "Depth of Ocean");
            ncvar.putAtt("units", zsgloO_units);
            ncvar.putAtt("sources", zsgloO_sources);
        auto ncvar2(ncio_blitz(ncio, mergemaskO, "MERGEMASK", "short", dims));
            ncvar.putAtt("description", "Identifies where merging changed the background fields");

    }

    if (errors.size() > 0) return -1;
    return 0;
}
