/** One-way coupler driver */

#include <icebin/modele/GCMCoupler_ModelE.hpp>

using namespace icebin;
using namespace ibmisc;


# Name of constants file written by ModelE
std::string const constants_fname = 'log/constants.nc'
std::string const constants_fname = 'log/constants.nc'

// =============================================================================
struct GCMOutputBundles {
    ibmisc::ArrayBundle<double,3> E;
};

class GCMOutputs {
public:
    GCMOutputBundles bundles;
    blitz::Array<double,3> &runo;
    blitz::Array<double,3> &eruno;
    blitz::Array<double,3> &deltah;
    blitz::Array<double,3> &massxfer;
    blitz::Array<double,3> &enthxfer;
    blitz::Array<double,3> &volxfer;
    blitz::Array<double,3> &gcm_bottom_senth;

    GCMOutputs(GCMOutputBundles &&_bundleE);
};

GCMOutputs::GCMOutputs(GCMOutputBUndles &&_bundles) :
    bundles(std::move(_bundles)),
    runo(bundles.E.array("runo")),
    eruno(bundles.E.array("eruno")),
    deltah(bundles.E.array("deltah")),
    massxfer(bundles.E.array("massxfer")),
    enthxfer(bundles.E.array("enthxfer")),
    volxfer(bundles.E.array("volxfer")),
    gcm_bottom_senth(bundles.E.array("gcm_bottom_senth"))
{}

ibmisc::ArrayBundle<double,3> gcm_outputs_bundleE()
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

class GCMInputs {
public:
    ibmisc::ArrayBundle<double,3> bundleE;
    ibmisc::ArrayBundle<double,2> bundleA;

    blitz::Array<double,2> &elevA;
    blitz::Array<double,3> &elevE;
    blitz::Array<double,3> &ice_top_senth;
    blitz::Array<double,2> & basal_frictional_heating;
    blitz::Array<double,2> & strain_heating;
    blitz::Array<double,2> & geothermal_flux;
    blitz::Array<double,2> & upward_geothermal_flux;
    blitz::Array<double,2> & calving_mass;
    blitz::Array<double,2> & calving_enth;
    blitz::Array<double,2> & basal_runoff_mass;
    blitz::Array<double,2> & basal_runoff_enth;
    blitz::Array<double,2> & internal_advection_mass;
    blitz::Array<double,2> & internal_advection_enth;
    blitz::Array<double,2> & epsilon_mass;
    blitz::Array<double,2> & epsilon_enth;

    GCMInputs(ibmisc::ArrayBundle<double,3> &&_bundleE);
};

GCMInputs::GCMInputs(GCMInputBundles &&_bundles) :
    bundles(std::move(_bundles)),
    elevA(bundles.A.array("elevA")),
    elevE(bundles.E.array("elevE")),
    ice_top_senth(bundles.E.array("ice_top_senth")),
    basal_frictional_heating(bundles.A.array("basal_frictional_heating")),
    strain_heating(bundles.A.array("strain_heating")),
    geothermal_flux(bundles.A.array("geothermal_flux")),
    upward_geothermal_flux(bundles.A.array("upward_geothermal_flux")),
    calving.mass(bundles.A.array("calving.mass")),
    calving.enth(bundles.A.array("calving.enth")),
    basal_runoff.mass(bundles.A.array("basal_runoff.mass")),
    basal_runoff.enth(bundles.A.array("basal_runoff.enth")),
    internal_advection.mass(bundles.A.array("internal_advection.mass")),
    internal_advection.enth(bundles.A.array("internal_advection.enth")),
    epsilon.mass(bundles.A.array("epsilon.mass")),
    epsilon.enth(bundles.A.array("epsilon.enth"))
{}

std::tuple<ibmisc::ArrayBundle<double,2>, ibmisc::ArrayBundle<double,3>> gcm_inputs_bundles()
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

    return ret;
}

// =============================================================================
struct GlobalBundles {
    ibmisc::ArrayBundle<double,3> Ed;
    ibmisc::ArrayBundle<int,3> Ei;
    ibmisc::ArrayBundle<double,2> Ad;
};

class Globals {
public:
    GlobalBundles bundles;

    blitz::Array<double, 3> &fhc;
    blitz::Array<int, 3> &underice;
    blitz::Array<double, 3> &elevE;
    blitz::Array<double, 2> &focean;
    blitz::Array<double, 2> &flake;
    blitz::Array<double, 2> &fgrnd;
    blitz::Array<double, 2> &fgice;
    blitz::Array<double, 2> &zatmo;

    Globals(GlobalBundles &&_bundles);
};

Globals::Globals() :
    bundles(std::move(_bundles)),
    fhc(bundles.Ed.array("fhc")),
    underice(bundles.Ei.array("underice")),
    elevE(bundles.Ed.array("elevE")),
    focean(bundles.Ad.array("focean")),
    flake(bundles.Ad.array("flake")),
    fgrnd(bundles.Ad.array("fgrnd")),
    fgice(bundles.Ad.array("fgice")),
    zatmo(bundles.Ad.array("zatmo"))
{}

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

    return ret;
}

// =============================================================================

struct Oneway {
    std::unique_ptr<GCMCoupler_ModelE> gcmce;

    int const im = 144;
    int const jm = 90;
    int nhc_gcm;

    Oneway::Oneway(int argc, char **argv)

    blitz::Array<double,3> add_gcm_inputE(
        std::string const &field, std::string const &units, std::string const &description, bool intitial);

    blitz::Array<double,2> add_gcm_inputA(
        std::string const &field, std::string const &units, std::string const &description, bool initial=false);

    blitz::Array<double,3> add_gcm_outputE(
        std::string const &field, std::string const &units, std::string const &description);


};

blitz::Array<double,3> Oneway::add_gcm_inputsE(ibmisc::ArrayBundle<double,3> &bundleE)
{
    blitz::Array<double,3> var_f(im,jm,nhc_gcm, blitz::fortranArray);
    gcmce_add_gcm_inpute(gcmce, var_f,
        field.c_str(), field.size(),
        units.c_str(), units.size(),
        description.c_str(), description.size());
    return var_f;
}

blitz::Array<double,2> Oneway::add_gcm_inputA(
    std::string const &field, std::string const &units, std::string const &description, bool initial)
{
    blitz::Array<double,2> var_f(im,jm, blitz::fortranArray);
    gcmce_add_gcm_inputa(gcmce, var_f,
        field.c_str(), field.size(),
        units.c_str(), units.size(),
        description.c_str(), description.size());
    return var_f;
}



blitz::Array<double,3> Oneway::add_gcm_outputE(
    std::string const &field, std::string const &units, std::string const &description)
{
    blitz::Array<double,3> var_f(im,jm,nhc_gcm, blitz::fortranArray);
    gcmce_add_gcm_outpute(gcmce, var_f,
        field.c_str(), field.size(),
        units.c_str(), units.size(),
        description.c_str(), description.size());
    return var_f;
}


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
            0);
        :
        gcmce_new(
            ModelEParams(),    // Dummy for now
            im, jm,
            1, 0, 1, 0,        // No part of the domain is in non-root
            MPI_Comm_c2f(comm),
            0);
    );

    // Figure out how many elevation classes we need.
    int nhc_gcm;
    int icebin_base_hc;
    int nhc_ice;
    gcmce_hc_params(gcmce, nhc_gcm, icebin_base_hc, nhc_ice);

    // Read constants from NetCDF file written during ModelE run
    // (in lieu of calling gcmce_set_constant())
    {NcIO ncio(consants_fname, 'r');
        gcmce->gcm_constants.read_nc(ncio.nc, "");
    }

    // -----------------------------------------------------------
    // gcmce_add_gcm_[input/output][a/e]()
    //
    // Allocate arrays used to communicate with coupler
    // And register them with IceBin Coupler
    GCMOutputs outputs(gcm_outputs_bundleE());
    outputs.bundles.E.allocate(
        blitz::shape(im,jm,nhc_gcm),
        {"im","jm","nhc_gcm"},
        blitz::fortranArray);
    for (size_t i=0; i<outputs.bundleE.index.size(); ++i) {
        auto ArrayBundle<double,3>::Meta const &meta(outputs.bundleE.data[i]);
        std::map<std::string, std::string> attr(meta.make_attr_map());

        std::string const &name(meta.name);
        std::string const &units(attr.at("units"));
        std::string const &description(attr.at("description"));
        gcmce_add_inpute(gcmce, meta.arr,
            name.c_str(), name.size(),
            units.c_str(), units.size(),
            description.c_str(), description.size());
    }

    GCMInputs inputs(gcm_inputs_bundles());
    inputs.bundles.A.allocate(
        blitz::shape(im,jm),
        {"im","jm"},
        blitz::fortranArray);
    for (size_t i=0; i<inputs.bundleA.index.size(); ++i) {
        auto ArrayBundle<double,2>::Meta const &meta(inputs.bundleA.data[i]);
        std::map<std::string, std::string> attr(meta.make_attr_map());

        std::string const &name(meta.name);
        std::string const &units(attr.at("units"));
        std::string const &description(attr.at("description"));
        auto initial_ii(attr.find("initial"));
        bool initial = (initial_ii == attr.end() ? true : atoi(initial_ii->second.c_str()));

        gcmce_add_inputa(gcmce, meta.arr,
            name.c_str(), name.size(),
            units.c_str(), units.size(),
            initial,
            description.c_str(), description.size());
    }

    inputs.bundles.E.allocate(
        blitz::shape(im,jm, nhc_gcm),
        {"im","jm","nhc_gcm"},
        blitz::fortranArray);
    for (size_t i=0; i<inputs.bundleE.index.size(); ++i) {
        auto ArrayBundle<double,3>::Meta const &meta(inputs.bundleE.data[i]);
        std::map<std::string, std::string> attr(meta.make_attr_map());

        std::string const &name(meta.name);
        std::string const &units(attr.at("units"));
        std::string const &description(attr.at("description"));
        auto initial_ii(attr.find("initial"));
        bool initial = (initial_ii == attr.end() ? true : atoi(initial_ii->second.c_str()));

        gcmce_add_inpute(gcmce, meta.arr,
            name.c_str(), name.size(),
            units.c_str(), units.size(),
            initial,
            description.c_str(), description.size());
    }

    // ----------------------------------------------------------------
    // gcmce_reference_globals()

    Globals 




    




    // Set up a single domain
    GCMCoupler_ModelE *gcmc

}
