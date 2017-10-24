#pragma once

#include <cstdlib>

#include <ibmisc/datetime.hpp>
#include <ibmisc/VarTransformer.hpp>
#include <ibmisc/DynArray.hpp>
#include <ibmisc/udunits2.hpp>
#include <ibmisc/ConstantSet.hpp>

#include <icebin/GCMParams.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/VarSet.hpp>
#include <icebin/multivec.hpp>

namespace ibmisc {
    class NcIO;
}

namespace icebin {

class GCMCoupler;
class GCMInput;    // formerly GCMCouplerOutput
class IceWriter;

class IceCoupler {
    friend class IceWriter;

public:
    BOOST_ENUM_VALUES( Type, int,
        (DISMAL)        (0)     // Demo Ice Sheet Model and LandIce
        (PISM)          (1)
        (ISSM)          (2)
        (WRITER)        (3)
    );
    const IceCoupler::Type type;

    /** Ordered specification of the variables (w/ units)
    to be passed IceBin->IceCoupler and IceCoupler->IceBin */
    enum IO {INPUT, OUTPUT, count};

public:
    GCMCoupler const *gcm_coupler;      // parent back-pointer
    IceRegridder const *ice_regridder;   // Set from gcm_coupler.
    std::string _name;              // Name of the ice sheet (in case ice_regridder is NULL)

    /** Place where we can write stuff related to this ice sheet */
    std::string output_dir;

    // Densified regridding matrix, and dimension, from previous call
    // Used to interpret GCM output
    std::unique_ptr<EigenSparseMatrixT> IvE0;
    // DenseArrayT<1> wIvE0;    // Provides the mask for I; for debugging only.
    SparseSetT dimE0;

    // Output of ice model from the last time we coupled.
    // Some of these values are needed for computation of ice_ivalsI
    // on the next coupling timestep.
    // TODO: Only keep variables we need here, instead of all of ice_ovalsI
    blitz::Array<double,2> ice_ovalsI;

    // [INPUT|OUTPUT] variables
    // List of fields this dynamic ice model takes for input / output.
    std::array<VarSet, IO::count> contract;

    // [INPUT|OUTPUT] variables
    // Defines IceBin-standard names for certain fiels in the contracts.
    // Allows IceBin to access these fields
    std::array<std::map<std::string,int>, IO::count> standard_names;

    // Linear combination transforming variables from:
    //      INPUT: gcm_output --> ice_input
    //     OUTPUT: ice_output --> gcm_input
    // (eg: T_ice = T_gcm + 273.15)
    ibmisc::VarTransformer var_trans_inE;
    std::array<ibmisc::VarTransformer, GridAE::count> var_trans_outAE;

    // Writers called to record the input and output seen by this IceCoupler
    std::array<std::unique_ptr<IceWriter>, 2> writer;

    // Current ice sheet elevation
    blitz::Array<double,1> elevI;
public:
    std::string const &name() const { return _name; }
    Grid const *gridI() { return ice_regridder->gridI.get(); }
    long nI() const { return ice_regridder->gridI->ndata(); }

    // ======================================================

    virtual ~IceCoupler();

public:
    // Lifecycle

    // ========= Called by
    //     **** PART 1: Allocate
    //     ATM_DRV.f: alloc_drv_atm()
    //     LISnow%allocate()
    //     GCMCoupler_ModelE: gcmce_new()
    //     GCMCoupler::ncread()
    //     IceCoupler.cpp: new_ice_coupler()
    IceCoupler(IceCoupler::Type _type);

    /** (1) Initialize any grid information, etc. from the IceSheet struct.
    @param vname_base Construct variable name from this, out of which to pull parameters from netCDF
    @param Opened handle on the IceBin config file (not IceBin grid file). */
    virtual void ncread(ibmisc::NcIO &ncio_config, std::string const &vname_sheet) {}


    // ========= Called by
    //     LANDICE_DRV.f: init_LI(istart_fixup)
    //     lisnow%cold_start()  (if cold start)
    //     lisheet%cold_start()  [currently missing...?]
    //     GCMCoupler::cold_start()
    //         <this>
    //     [calls IceCoupler::cold_start()]
    /** (2) Event handler to let IceCouplers know the start time is (finally) set.
    This method must call where appropriate:
        contracts::setup(*gcm_coupler, *this);
    */
    void cold_start(
        ibmisc::Datetime const &time_base,
        double time_start_s);

    /** Called by cold_start() */
    virtual void _cold_start(
        ibmisc::Datetime const &time_base,
        double time_start_s) = 0;

    void print_contracts();

    // ========= Called by

protected:

    /** Called in just one place.  Transforms GCM output to Ice Model
    input using general means.  This function is customized by
    setting reconstruct_ice_ivalsI below. */
    blitz::Array<double,2> construct_ice_ivalsI(
        blitz::Array<double,2> const &gcm_ovalsE0,
        std::vector<std::pair<std::string, double>> const &scalars,
        double dt,
        ibmisc::TmpAlloc &tmp);

public:
    /** A "virtual function" used to customize construct_ice_ivalsI().
    This defaults to NOP, and is set by the coupling contract. */
    std::function<void(blitz::Array<double,2> &, double)> reconstruct_ice_ivalsI;

    /** (4) Run the ice model for one coupling timestep.
    @param time_s Seconds since GCMParams::time_base.  Helps with debugging.
    @param index Index of each input grid value in ivalsI.
    @param ivalsI The values themselves (sparse representation).
           Their meaning (SMB, T, etc) is determined
           by the place in the array, as specified by the appropriate
           INPUT contract for this ice model.
    @param am_i_root
        Call with true if calling from MPI root; false otherwise.
        The core coupling/regridding computation only runs on root.
        But other MPI ranks need to go along for the ride, assuming that
        the ice model uses MPI. */
    void couple(
        double time_s,
        // Values from GCM, passed GCM -> Ice
        VectorMultivec const &gcm_ovalsE,
        GCMInput &out,    // Accumulate matrices here...
        bool do_run);

    /** (4.1) @param index Index of each grid value.
    @param time_s Time since start of simulation, in seconds
    @param do_run True if we are to actually run (otherwise just return ice_ovalsI from current state) */
    virtual void run_timestep(double time_s,
        blitz::Array<double,2> const &ice_ivalsI,
        blitz::Array<double,2> &ice_ovalsI,
        bool run_ice) = 0;

};      // class IceCoupler
// =========================================================

extern
std::unique_ptr<IceCoupler> new_ice_coupler(ibmisc::NcIO &ncio,
    std::string const &vname, std::string const &sheet_name,
    GCMCoupler const *_gcm_coupler);


// =========================================================
class IceWriter
{
    IceCoupler const *ice_coupler;

    /** Description of the fields we're writing */
    VarSet const *contract;

    // The output file we are writing to...
    std::string fname;

    // Used for lazy opening of output file
    bool file_initialized = false;

    // Dimensions to use when writing to netCDF
    std::vector<std::string> dim_names;
    std::vector<size_t> cur;        // Base index to write in netCDF
    std::vector<size_t> counts;

public:
    IceWriter(
        IceCoupler const *_ice_coupler,
        VarSet const *_contract,
        std::string const &_fname);

    void write(double time_s,
        blitz::Array<double,2> const &valsI);    // valsI[nI, nvars]

private:
    void init_file();

};

}
