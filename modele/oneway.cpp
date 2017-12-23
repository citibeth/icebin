/** One-way coupler driver */

#include <icebin/modele/GCMCoupler_ModelE.hpp>
#include <ibmisc/bundle.hpp>
#include <everytrace.h>

using namespace icebin;
using namespace icebin::modele;
using namespace ibmisc;

// Name of constants file written by ModelE
std::string const constants_fname = "log/constants.nc";
std::string const TOPO_fname = "TOPO";
std::string const forcing_fname = "gcm-out-19500305.nc";

// =============================================================================
struct GCMOutputBundles {
    ibmisc::ArrayBundle<double,3> E;
};

GCMOutputBundles gcm_outputs_bundles()
{
    GCMOutputBundles bundles;

    bundles.E.add("runo", {
        "units", "kg m-2 s-1",
        "description", "Downward water flux through bottom layer"
    });
    bundles.E.add("eruno", {
        "units", "W m-2",
        "description", "Enthalpy of downward water flux through bottom layer"
    });
    // auto trruno(add_gcm_outputE(this, this%snogli%trruno, "trruno"));
    bundles.E.add("deltah", {
        "units", "W m-2",
        "description", "Enthalpy change of 'borrowed' layer"
    });
    bundles.E.add("massxfer", {
        "units", "kg m-2 s-1",
        "description", "Mass of ice being transferred Stieglitz --> Icebin"
    });
    bundles.E.add("enthxfer", {
        "units", "W m-2",
        "description", "Enthlpy of ice being transferred Stieglitz --> Icebin"
    });
    // auto trxfer(add_gcm_outputE(this, this%snogli%trxfer, "trxfer"));
    bundles.E.add("volxfer", {
        "units", "m^3 m-2 s-1",
        "description", "Volume of ice being transferred Stieglitz --> Icebin"
    });

    bundles.E.add("gcm_bottom_senth", {
        "units", "J kg-1",
        "description", "Sepcific Enthalpy at bottom of GCM's snow/firn model"
    });

    return bundles;
}
// =============================================================================
struct GCMInputBundles {
    ibmisc::ArrayBundle<double,2> A;
    ibmisc::ArrayBundle<double,3> E;
};

GCMInputBundles gcm_inputs_bundles()
{
    GCMInputBundles bundles;

    // --------- Dynamic Ice Model State
    bundles.E.add("elevA", {
        "units", "m",
        "description", "ice upper surface elevation",
        "initial", "1"
    });
    bundles.E.add("elevE", {
        "units", "m",
        "description", "ice upper surface elevation",
        "initial", "1"
    });
    bundles.E.add("ice_top_senth", {
        "units", "J kg-1",
        "description", "Specific enthalpy at top of ice model's ice sheet",
        "initial", "1"
    });

    // ---------- Heat Flux Ice Model Outputs
    bundles.A.add("basal_frictional_heating", {
        "units", "W m-2",
        "description", "Frictional heating at base of ice sheet"
    });
    bundles.A.add("strain_heating", {
        "units", "W m-2",
        "description", "Heating from internal friciton"
    });
    bundles.A.add("geothermal_flux", {
        "units", "W m-2",
        "description", "Heat flow between ice sheet and solid earth. ???"
    });
    bundles.A.add("upward_geothermal_flux", {
        "units", "W m-2",
        "description", "Heat flow between ice sheet and solid earth. ???"
    });


    // ----------- Mass Transfer Flux Outputs
    bundles.A.add("calving.mass", {
        "units", "kg m-2 s-1",
        "description", "Calving rate for grid cells containing a calving front."
    });
    bundles.A.add("calving.enth", {
        "units", "W m-2",
        "description", "Calving rate for grid cells containing a calving front."
    });
    bundles.A.add("basal_runoff.mass", {
        "units", "kg m-2 s-1",
        "description", "Basal melting of grounded ice"
    });
    bundles.A.add("basal_runoff.enth", {
        "units", "W m-2",
        "description", "Basal melting of grounded ice"
    });
    bundles.A.add("internal_advection.mass", {
        "units", "kg m-2 s-1",
        "description", "Horizontal advection due to ice dynamics"
    });
    bundles.A.add("internal_advection.enth", {
        "units", "W m-2",
        "description", "Horizontal advection due to ice dynamics"
    });
    bundles.A.add("epsilon.mass", {
        "units", "kg m-2 s-1",
        "description", "Changes not otherwise accounted for"
    });
    bundles.A.add("epsilon.enth", {
        "units", "W m-2",
        "description", "Changes not otherwise accounted for"
    });

    return bundles;
}

// =============================================================================
struct GlobalBundles {
    ibmisc::ArrayBundle<double,3> Ed;
    ibmisc::ArrayBundle<int,3> Ei;
    ibmisc::ArrayBundle<double,2> Ad;
};

GlobalBundles global_bundles()
{
    GlobalBundles bundles;

    bundles.Ed.add("fhc", {
        "units", "1",
        "description", "Contribution of elevation class to overall Atmosphere grid cell."
    });
    bundles.Ei.add("underice", {
        "units", "1=ICEBIN,2=NOTHING",
        "description", "Indicates where a local ice model vs. push-down stead-state resides."
    });
    bundles.Ed.add("elevE", {
        "units", "m",
        "description", "Elevation of ice surface"
    });
    bundles.Ad.add("focean", {
        "units", "1",
        "description", "Fraction of Atmosphere grid cell that is ocean."
    });
    bundles.Ad.add("flake", {
        "units", "1",
        "description", "Fraction of Atmosphere grid cell that is ocean."
    });
    bundles.Ad.add("focean", {
        "units", "1",
        "description", "Fraction of Atmosphere grid cell that is lake."
    });
    bundles.Ad.add("fgrnd", {
        "units", "1",
        "description", "Fraction of Atmosphere grid cell that is ice-free ground."
    });
    bundles.Ad.add("fgice", {
        "units", "1",
        "description", "Fraction of Atmosphere grid cell that is land ice."
    });
    bundles.Ad.add("zatmo", {
        "units", "m",
        "description", "Elevation of bottom of the atmosphere."
    });

    return bundles;
}

// =============================================================================

struct Oneway {
    std::unique_ptr<GCMCoupler_ModelE> gcmce;

    int const im = 144;
    int const jm = 90;
    int nhc_gcm;

    Oneway(int argc, char **argv);

};



Oneway::Oneway(int argc, char **argv)
{
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the communicator
    MPI_Comm comm = MPI_COMM_WORLD;

    // Get the number of processes
    int world_size;
    MPI_Comm_size(comm, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
 
    // Print a hello world message
    printf("Hello world from processor %s, rank %d"
        " out of %d processors\n",
        processor_name, world_rank, world_size);
 
    // -----------------------------------

    // Initialize.  Put the entire domain in root, to make our driver job easier.
    // NOTE: gcmce_new() reads from ./config/icebin.nc
    gcmce.reset(world_rank == 0 ?
        gcmce_new(
            ModelEParams(),    // Dummy for now
            im, jm,
            1, im, 1, jm,        // The entire domain is in root
            MPI_Comm_c2f(comm),
            0)
        :
        gcmce_new(
            ModelEParams(),    // Dummy for now
            im, jm,
            1, 0, 1, 0,        // No part of the domain is in non-root
            MPI_Comm_c2f(comm),
            0)
    );

    // Figure out how many elevation classes we need.
    int nhc_gcm;
    int icebin_base_hc;
    int nhc_ice;
    gcmce_hc_params(&*gcmce, nhc_gcm, icebin_base_hc, nhc_ice);

    // Read constants from NetCDF file written during ModelE run
    // (in lieu of calling gcmce_set_constant())
    {NcIO ncio(constants_fname, 'r');
        gcmce->gcm_constants.read_nc(ncio.nc, "");
    }

    // -----------------------------------------------------------
    // gcmce_add_gcm_[input/output][a/e]()
    //
    // Allocate arrays used to communicate with coupler
    // And register them with IceBin Coupler
    GCMOutputBundles outputs(gcm_outputs_bundles());
    outputs.E.allocate(
        {im,jm,nhc_gcm},
        {"im","jm","nhc_gcm"},
        true, blitz::fortranArray);
    for (size_t i=0; i<outputs.E.index.size(); ++i) {
        ArrayBundle<double,3>::Data &meta(outputs.E.data[i]);
//        std::map<std::string, std::string> attr(meta.meta.attr.begin(), meta.meta.attr.end());
        std::map<std::string, std::string> attr(meta.meta.make_attr_map());

        std::string const &name(meta.meta.name);
        std::string const &units(attr.at("units"));
        std::string const &description(attr.at("description"));
        auto initial_ii(attr.find("initial"));
        bool initial = (initial_ii == attr.end() ? true : atoi(initial_ii->second.c_str()));

        auto arr_f(f90array(meta.arr));
        gcmce_add_gcm_inpute(&*gcmce, arr_f,
            name.c_str(), name.size(),
            units.c_str(), units.size(),
            initial,
            description.c_str(), description.size());
    }

    GCMInputBundles inputs(gcm_inputs_bundles());
    inputs.A.allocate(
        {im, jm},
        {"im","jm"},
        true, blitz::fortranArray);
    for (size_t i=0; i<inputs.A.index.size(); ++i) {
        ArrayBundle<double,2>::Data &meta(inputs.A.data[i]);
        std::map<std::string, std::string> attr(meta.meta.make_attr_map());

        std::string const &name(meta.meta.name);
        std::string const &units(attr.at("units"));
        std::string const &description(attr.at("description"));
        auto initial_ii(attr.find("initial"));
        bool initial = (initial_ii == attr.end() ? true : atoi(initial_ii->second.c_str()));

        auto arr_f(f90array(meta.arr));
        gcmce_add_gcm_inputa(&*gcmce, arr_f,
            name.c_str(), name.size(),
            units.c_str(), units.size(),
            initial,
            description.c_str(), description.size());
    }

    inputs.E.allocate(
        {im,jm, nhc_gcm},
        {"im","jm","nhc_gcm"},
        true, blitz::fortranArray);
    for (size_t i=0; i<inputs.E.index.size(); ++i) {
        ArrayBundle<double,3>::Data &meta(inputs.E.data[i]);
        std::map<std::string, std::string> attr(meta.meta.make_attr_map());

        std::string const &name(meta.meta.name);
        std::string const &units(attr.at("units"));
        std::string const &description(attr.at("description"));
        auto initial_ii(attr.find("initial"));
        bool initial = (initial_ii == attr.end() ? true : atoi(initial_ii->second.c_str()));

        auto arr_f(f90array(meta.arr));
        gcmce_add_gcm_inpute(&*gcmce, arr_f,
            name.c_str(), name.size(),
            units.c_str(), units.size(),
            initial,
            description.c_str(), description.size());
    }

    // ----------------------------------------------------------------
    // gcmce_reference_globals()

    GlobalBundles globals;
    globals.Ed.allocate(
        {im,jm, nhc_gcm},
        {"im","jm","nhc_gcm"},
        true,
        blitz::fortranArray);
    globals.Ei.allocate(
        {im,jm, nhc_gcm},
        {"im","jm","nhc_gcm"},
        true,
        blitz::fortranArray);
    globals.Ad.allocate(
        {im,jm},
        {"im","jm"},
        true,
        blitz::fortranArray);

    gcmce_reference_globals(&*gcmce,
        f90array(globals.Ed.array("fhc")),
        f90array(globals.Ei.array("underice")),
        f90array(globals.Ed.array("elevE")),
        f90array(globals.Ad.array("focean")),
        f90array(globals.Ad.array("flake")),
        f90array(globals.Ad.array("fgrnd")),
        f90array(globals.Ad.array("fgice")),
        f90array(globals.Ad.array("zatmo"))
    );

    // -----------------------------------------------------------------
    // Read TOPO file into GlobalBundles
    {NcIO ncio(TOPO_fname, 'r');
        globals.Ed.ncio(ncio, {}, "", "double");
        globals.Ei.ncio(ncio, {}, "", "int");
        globals.Ad.ncio(ncio, {}, "", "double");
    }


    // -----------------------------------------------------------------
    // Read a single forcing to use over and over again...

    int itime=0;
    {NcIO ncio(forcing_fname, 'r');
        outputs.E.ncio_partial(
            ncio, {}, "", "double",
            {}, {itime, 0, 0, 0}, {3,2,1});
    }


    // Set up a single domain
    GCMCoupler_ModelE *gcmc;

}


int main(int argc, char **argv)
{
    everytrace_init();
    Oneway ow(argc, argv);
}

