/*
INPUT:

GCM grid spec  (string: grid name)

Ice (hi-res) grid spec
FGICE, elev --> elevmaskI
    --> AbbrGrid, only of ice-covered grid cells

elevations for ECs


OUTPUT:

unscaled matrices IvA, IvE, AvE


Run with:
    ulimit -v 8000000

*/


#include <string>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <tclap/CmdLine.h>

#include <prettyprint.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/stdio.hpp>
#include <ibmisc/filesystem.hpp>
#include <spsparse/SparseSet.hpp>

//#include <icebin/Grid.hpp>
//#include <icebin/AbbrGrid.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/modele/GCMRegridder_ModelE.hpp>
#include <icebin/modele/grids.hpp>
#include <icebin/modele/hntr.hpp>

#include <icebin/gridgen/GridGen_LonLat.hpp>

using namespace std;
using namespace ibmisc;
using namespace icebin;
using namespace netCDF;
using namespace spsparse;

// ==========================================================
struct ParseArgs {
    HntrSpec hspecO;    // Name of Hntr Spec for Ocean grid
    HntrSpec hspecI;    // Name of Hntr Spec for Ice Grid

    std::string nc_fname;
    std::string fgiceI_vname;
    std::string elevI_vname;

//    std::string nc foceanI_fname;
    std::string topoO_fname;
//    std::string foceanI_vname;
//    std::string 

    std::string ofname;
    std::array<double,2> ec_range;    // Lowest and highest elevation classes [m]
    double ec_skip;                    // Distance between elevation classes [m]
    bool scale;
    bool const correctA = false;    // Only needed with projected I grids (and then not really)
    std::array<double,3> sigma;

    double eq_rad;        // Radius of earth; see ModelE code    

    ParseArgs(int argc, char **argv);
};


ostream& operator<<(ostream& os, ParseArgs const &args)
{  
    return os << "ParseArgs(" << endl
        << "    hspecO: " << args.hspecO.im << "x" << args.hspecO.jm << endl
        << "    hspecI: " << args.hspecI.im << "x" << args.hspecI.jm << endl
        << "    nc_fname: " << args.nc_fname << " -- " << args.fgiceI_vname << " -- " << args.elevI_vname << endl
        << "    topoO_fname: " << args.topoO_fname << endl
        << "    ofname: " << args.ofname << endl
        << "    ec_range: " << args.ec_range << "  ec_skip=" << args.ec_skip << endl
        << "    scale: " << (args.scale ? "true" : "false") << endl
        << "    sigma: " << args.sigma;
}  

template<class DestT>
static std::vector<DestT> parse_csv(std::string scsv_str)
{
    // Parse to vector of strings
    std::vector<std::string> scsv;
    boost::algorithm::split(scsv,  scsv_str, boost::is_any_of(","));


    std::vector<DestT> ret;
    for (std::string &s : scsv) {
        stringstream myString(s);
        DestT val;
        myString >> val;
        ret.push_back(val);
    }


    return ret;
}

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

        TCLAP::UnlabeledValueArg<std::string> shspecO_a(
            "gridO", "Name of Ocean grid (eg: g1qx1)",
            true, "g1qx1", "atm. grid name", cmd);

        TCLAP::UnlabeledValueArg<std::string> shspecI_a(
            "gridI", "Name of Ice grid (eg: g1mx1m)",
            true, "g1mx1m", "ice grid name", cmd);

        TCLAP::UnlabeledValueArg<std::string> nc_fname_a(
            "elevmaskI-fname",
                "NetCDF file containing ice mask and elevation (1 where there is ice)",
            true, "etopo1_ice_g1m.nc", "mask filename", cmd);

        TCLAP::UnlabeledValueArg<std::string> topoO_fname_a(
            "topoO-fname",
                "ModelE TOPO file on the Ocean grid.  Need FOCEAN and FOCEANF",
            true, "topoo.nc", "focean filename", cmd);





        TCLAP::ValueArg<std::string> foceanI_vname_a("n", "focean",
            "Name of NetCDF variable containing ice focean (1 where there is ice)",
            false, "FGICE1m", "focean var name", cmd);

        TCLAP::ValueArg<std::string> fgiceI_vname_a("m", "mask",
            "Name of NetCDF variable containing ice mask (1 where there is ice)",
            false, "FGICE1m", "mask var name", cmd);

        TCLAP::ValueArg<std::string> elevI_vname_a("e", "elev",
            "Name of NetCDF variable containing elevation [m]",
            false, "ZICETOP1m", "elevation var name", cmd);

        TCLAP::ValueArg<std::string> ec_a("E", "elev-classes",
            "Elevations [m] for the elevation classes: lowest,highest,skip",
            false, "-100,3700,200", "elevations", cmd);

        TCLAP::ValueArg<std::string> ofname_a("o", "output",
            "Output filename (NetCDF) for ECs",
            false, "global_ec.nc", "mask var name", cmd);

        TCLAP::SwitchArg raw_a("r", "raw",
             "Produce raw (unscaled) matrices?",
             cmd, false);

        TCLAP::ValueArg<std::string> sigma_a("g", "sigma",
            "Scaling factor: x,y,z",
            false, "0,0,0", "scaling factor", cmd);

        TCLAP::ValueArg<double> eq_rad_a("R", "radius",
            "Radius of the earth",
            false, modele::EQ_RAD, "earth radius", cmd);


        // Not needed for spherical grids
        // TCLAP::SwitchArg correctA_a("c", "correct",
        //      "Correct for area changes due to projection?",
        //      cmd, false);

        // Parse the argv array.
        cmd.parse( argc, argv );

        // Get the value parsed by each arg.
        hspecO = *modele::grids.at(shspecO_a.getValue());
        hspecI = *modele::grids.at(shspecI_a.getValue());
        nc_fname = nc_fname_a.getValue();
        fgiceI_vname = fgiceI_vname_a.getValue();
        elevI_vname = elevI_vname_a.getValue();
        topoO_fname = topoO_fname_a.getValue();

        // Parse elevation classes...
        auto _ec(parse_csv<double>(ec_a.getValue()));
        if (_ec.size() < 2 || _ec.size() > 3) (*icebin_error)(-1,
            "--ec '%s' must have just two or three values");
        ec_range[0] = _ec[0];
        ec_range[1] = _ec[1];
        ec_skip = (_ec.size() == 3 ? _ec[2] : 1);

        ofname = ofname_a.getValue();
        scale = !raw_a.getValue();

        eq_rad = eq_rad_a.getValue();

        // Parse sigma
        auto _sigma(parse_csv<double>(sigma_a.getValue()));
        if (sigma.size() != 3) (*icebin_error)(-1,
            "sigma must have exactly three elements");
        for (int i=0; i<_sigma.size(); ++i) sigma[i] = _sigma[i];

    } catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        exit(1);
    }
}

// ==========================================================

class ExchAccum {
    ExchangeGrid &exgrid;
    blitz::Array<double,1> const &elevmaskI;
    SparseSet<long,int> &dimO;        // Dimension is created but not used
    SparseSet<long,int> &dimI;        // Dimension is created but not used
public:
    ExchAccum(
        ExchangeGrid &_exgrid,
        blitz::Array<double,1> const &_elevmaskI,
        SparseSet<long,int> &_dimO,
        SparseSet<long,int> &_dimI)
    : exgrid(_exgrid), elevmaskI(_elevmaskI), dimO(_dimO), dimI(_dimI) {}

    void add(std::array<int,2> const &index, double area)
    {
        auto const iO(index[0]);
        auto const iI(index[1]);
        if (!std::isnan(elevmaskI(iI))) {
            // Save as sparse indexing, as required by IceRegridder::init()
            exgrid.add(index, area);
            dimO.add_dense(iO);
            dimI.add_dense(iI);
        }
    }
};

static double const NaN = std::numeric_limits<double>::quiet_NaN();


int main(int argc, char **argv)
{
    everytrace_init();
    ParseArgs args(argc, argv);
    std::cout << args << endl;
    EnvSearchPath files("MODELE_FILE_PATH");

    ExchangeGrid aexgrid;    // Put our answer in here

    auto const &hspecO(args.hspecO);
    auto const &hspecI(args.hspecI);
    modele::Hntr hntr(17.17, hspecO, hspecI);

#if 0
    // Load the fractional ocean mask (based purely on ice extent)
    blitz::Array<double,2> foceanO(hspecO.jm, hspecO.im);    // called FOCEANF in make_topoo
    {
        blitz::Array<int16_t,2> foceanI(hspecI.jm, hspecI.im);
        {NcIO ncio(files.locate(args.focean_fname), 'r');
            ncio_blitz(ncio, foceanI, args.foceanI_vname, "short", {});
        }


        // Regrid O <- I
        blitz::Array<double, 2> wt1m(const_array(shape(hspecO.jm, hspecO.im), 1.0));
        hntr.regrid<double,int16_t,double>(wt1m, foceanI, foceanO, true);
    }
#endif


    // Load the ice mask
    blitz::Array<double,2> elevmaskI(hspecI.jm, hspecI.im);
    {
        auto fname(files.locate(args.nc_fname));
        printf("---- Reading elevmaskI: %s\n", fname.c_str());

        // Read in ice extent and elevation
        blitz::Array<int16_t,2> fgiceI(hspecI.jm, hspecI.im);    // 0 or 1
        blitz::Array<int16_t,2> elevI(hspecI.jm, hspecI.im);
        {NcIO ncio(fname, 'r');
            ncio_blitz(ncio, fgiceI, args.fgiceI_vname, "short", {});
            ncio_blitz(ncio, elevI, args.elevI_vname, "short", {});
        }

        // Combine into IceBin-standard elevmaskI
int nice=0;
        for (int j=0; j<hspecI.jm; ++j) {
        for (int i=0; i<hspecI.im; ++i) {
            elevmaskI(j,i) = (fgiceI(j,i) ? elevI(j,i) : NaN);
if (fgiceI(j,i)) ++nice;
        }}
printf("nice = %d\n", nice);
    }

    // -------------------------------------------------------------
    printf("---- Computing overlaps\n");

    // Compute overlaps for cells with ice
    SparseSet<long,int> _dimO;    // Only include O grid cells with ice
    SparseSet<long,int> _dimI;    // Only include I grid cells with ice
    hntr.overlap(ExchAccum(aexgrid, reshape1(elevmaskI), _dimO, _dimI), args.eq_rad);
printf("|dimO|=%d\n", _dimO.dense_extent());
//return 0;
#if 0
    {NcIO ncio(args.ofname, 'w');
        aexgrid.ncio(ncio, "aexgrid");
    }
#endif

    // -------------------------------------------------------------
    printf("---- Creating gcmO\n");

    // Turn HntrSpec --> GridSpec
    GridSpec_LonLat specO(make_grid_spec(hspecO, false, 1, args.eq_rad));
    GridSpec_LonLat specI(make_grid_spec(hspecI, false, 1, args.eq_rad));

    // Realize O grid for relevant gridcells
    auto agridO(make_abbr_grid("Ocean", specO, std::move(_dimO)));
printf("nO = %ld\n", agridO.dim.sparse_extent());

    // Set up elevation classes    
    std::vector<double> hcdefs;
    for (double elev=args.ec_range[0]; elev <= args.ec_range[1]; elev += args.ec_range[2]) {
        hcdefs.push_back(elev);
    }

    // Create standard GCMRegridder for O <--> I
    std::unique_ptr<GCMRegridder_Standard> gcmO(new GCMRegridder_Standard);
    gcmO->init(
        std::move(agridO), std::move(hcdefs),
        Indexing({"O", "HC"}, {0,0}, {agridO.dim.sparse_extent(), hcdefs.size()}, {1,0}),
        args.correctA);

    // --------------------------------------------------
    // Create IceRegridder for I and add to gcmO
    auto ice(new_ice_regridder(IceRegridder::Type::L0));
    auto agridI(make_abbr_grid("Ice", specI, std::move(_dimI)));
    ice->init("globalI", gcmO->agridA, nullptr,
        std::move(agridI), std::move(aexgrid),
        InterpStyle::Z_INTERP);    // You can use different InterpStyle if you like.

    gcmO->add_sheet(std::move(ice));

    // --------------------------------------------------
    printf("---- Creating gcmA\n");

    // Create a mismatched regridder, to mediate between different ice
    // extent of GCM vs. IceBin
    modele::GCMRegridder_ModelE gcmA(std::shared_ptr<GCMRegridder>(gcmO.release()));

    // Load the fractional ocean mask (based purely on ice extent)
    {auto fname(files.locate(args.topoO_fname));

printf("nO = %ld (%d %d)\n", gcmA.gcmO->nA(), hspecO.jm, hspecO.im);

        blitz::Array<double,2> foceanO(hspecO.jm, hspecO.im);    // called FOCEAN in make_topoo
        blitz::Array<double,2> foceanfO(hspecO.jm, hspecO.im);    // called FOCEANF in make_topoo

        printf("---- Reading FOCEAN: %s\n", fname.c_str());
        NcIO ncio(fname, 'r');
        ncio_blitz(ncio, foceanO, "FOCEAN", "double", {});
        ncio_blitz(ncio, foceanfO, "FOCEANF", "double", {});


printf("FF1\n");
    gcmA.foceanAOp = reshape1(foceanfO);        // COPY
printf("FF2\n");
    gcmA.foceanAOm = reshape1(foceanO);         // COPY
printf("FF3\n");
    }
printf("FF4\n");

    RegridMatrices rm(gcmA.regrid_matrices(0, reshape1(elevmaskI)));
printf("FF5\n");
    elevmaskI.free();    // No longer needed, and it is BIG
printf("FF6\n");

    // ---------- Generate and store the matrices
    // Use the mismatched regridder to create desired matrices and save to file
    RegridMatrices::Params params(args.scale, args.correctA, args.sigma);
    SparseSet<long,int> dimA, dimI, dimE;

    {NcIO ncio(args.ofname, 'w');
        printf("---- Generating AvI\n");
        auto mat(rm.matrix("AvI", {&dimA, &dimI}, params));
printf("dimA: %ld %ld\n", dimA.dense_extent(), dimA.sparse_extent());
printf("dimI: %ld %ld\n", dimI.dense_extent(), dimI.sparse_extent());
        mat->ncio(ncio, "AvI", {"dimA", "dimI"});
    }

#if 0
    {NcIO ncio(args.ofname, 'a');
        printf("---- Generating EvI\n");
        auto mat(rm.matrix("EvI", {&dimE, &dimI}, params));
printf("dimE: %ld %ld\n", dimE.dense_extent(), dimE.sparse_extent());
printf("dimI: %ld %ld\n", dimI.dense_extent(), dimI.sparse_extent());
printf("mat size: %ld\n", (long)(mat->M->nonZeros()));
        mat->ncio(ncio, "EvI", {"dimE", "dimI"});
    }
#endif

#if 0
    {NcIO ncio(args.ofname, 'a');
        printf("---- Generating IbE\n");
        auto mat(rm.matrix("IvE", {&dimI, &dimE}, params));
        mat->ncio(ncio, "IvE", {"dimI", "dimE"});
    }
    {NcIO ncio(args.ofname, 'a');
        printf("---- Generating IvA\n");
        auto mat(rm.matrix("IvA", {&dimI, &dimA}, params));
        mat->ncio(ncio, "IvA", {"dimI", "dimA"});
    }

    {NcIO ncio(args.ofname, 'a');
        printf("---- Generating AvE\n");
        auto mat(rm.matrix("AvE", {&dimA, &dimE}, params));
        mat->ncio(ncio, "AvE", {"dimA", "dimE"});
    }
#endif

    // Store the dimensions
    printf("---- Storing Dimensions\n");
    {NcIO ncio(args.ofname, 'a');
        dimA.ncio(ncio, "dimA");
        dimE.ncio(ncio, "dimE");
        dimI.ncio(ncio, "dimI");
    }

    printf("Done!\n");
    return 0;
}
