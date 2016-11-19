#pragma once

#include <cstdlib>

#include <ibmisc/VarTransformer.hpp>
#include <ibmisc/DynArray.hpp>
#include <ibmisc/udunits2.hpp>
#include <ibmisc/ConstantSet.hpp>

#include <icebin/GCMParams.hpp>
#include <icebin/GCMPerIceSheetParams.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/VarSet.hpp>

namespace icebin {

class GCMCoupler;
class IceCoupler_Writer;

class IceCoupler {
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
    enum IO {INPUT, OUTPUT};

//  friend class GCMCoupler;
//    friend class IceCoupler_Writer;
//  friend std::unique_ptr<IceCoupler> new_ice_model(IceCoupler::Type type,
//      GCMCoupler const *_coupler, IceRegridder const *_sheet);

public:
    GCMCoupler const *coupler;      // parent back-pointer
    IceRegridder *regridder;   // This is not const; see IceCoupler::update_elevI()
    std::unique_ptr<WeightedSparse> IvE;    // Regridding matrix made from regridder

    std::array<VarSet, 2> contract;     // [INPUT|OUTPUT]
    std::array<ibmisc::VarTransformer, 2> var_transformer;

    // Parameters provided by the GCM, to inform the coupling
    std::unique_ptr<GCMPerIceSheetParams> gcm_per_ice_sheet_params;

    // Writers called to record the input and output seen by this IceCoupler
    std::unique_ptr<IceCoupler_Writer> _iwriter, _owriter;

public:
    std::string const &name() { return sheet->name(); }
    Grid const *gridI() { return &*sheet->gridI; }

    /** Constants obtained from the GCM */
    ibmisc::ConstantSet ice_constants;

    


    /** Placeholder for additional coupling contracts that had to be allocated. */
//  std::vector<std::unique_ptr<giss::CouplingContract>> _extra_contracts;

    // --------------------------------------------
    // Buffers used to receive ice model output, and regrid it.

    /** Direct output from the ice model (on the ice grid).
    This is allocated (for the ice model ROOT node only) inside
    run_timestep() */
    std::vector<blitz::Array<double,1>> ice_ovals_I;

    /** Input to the GCM, but on the ice grid */
    std::vector<blitz::Array<double,1>> gcm_ivals_I;
    // ======================================================

    /** @return the iwriter associated with this IceCoupler, if the
    IceCoupler is NOT responsible for calling the writer itself.
    If the IceCoupler WILL call the writer, returns NULL. */
    virtual IceCoupler_Writer *iwriter() { return _iwriter.get(); }
    virtual IceCoupler_Writer *owriter() { return _owriter.get(); }


    /** Tells whether we are running on the root node of the ice model. */
    virtual bool am_i_root() const;

    /** Allocate vectors in preparation of calling an ice model (ROOT only). */
    void allocate_ice_ovals_I();

    /** Allocate in preparation of var transformations (but not regridding yet) (ROOT only) */
    void allocate_gcm_ivals_I();

    /** Free portions not needed after finished calling ice model and
    applying variable transform.  This will be variables desired on
    anything other than the ELEVATION grid. (ROOT only) */
    void free_ice_ovals_I();

    /** Free all memory used by this.  Called when we're done with a coupling timestep. (ROOT only) */
    void free_ovals_ivals_I();

    /** Allocates and sets gcm_ivals_I variable */
    void set_gcm_inputs(unsigned int mask);

    // --------------------------------------------
#if 0
    /** Allocate a new giss::CouplingContract, with the same lifetime as this IceCoupler. */
    VarSet *new_VarSet();
#endif

    IceCoupler(IceCoupler::Type _type) : type(_type) {}
    virtual ~IceCoupler();

    long ndata() const { return sheet->gridI->ndata(); }

    // --------------------------------------------------

public:

    /** Initialize any grid information, etc. from the IceSheet struct.
    @param vname_base Construct variable name from this, out of which to pull parameters from netCDF */
    virtual void ncread(ibmisc::NcIO &ncio, std::string const &vname_sheet);

    /** Event handler to let IceCouplers know the start time is (finally) set */
    virtual void start_time_set() {}

    /** Sets elevI based on the latest state from the ice model. */
    virtual void _update_elevI() = 0;

    RegridMatrices update_elevI()
    {
        // Get new ice sheet elevations and mask
        _update_elevI();
        // Use that to re-do our regridding matrix.
        rm = RegridMatrices(regridder);
        auto IvE(re.regrid("IvE", true, true);
        IvE = 
        return rm;
    }

    /** Run the ice model for one coupling timestep.
    @param time_s Seconds since GCMParams::time_base.  Helps with debugging.
    @param index Index of each input grid value in ivalsI.
    @param ivalsI The values themselves (sparse representation).
           Their meaning (SMB, T, etc) is determined
           by the place in the array, as specified by the appropriate
           INPUT contract for this ice model.
    */
    virtual void run_timestep(double time_s,
        blitz::Array<int,1> const &indices,
        std::vector<blitz::Array<double,1>> const &ivalsI)
    {
        (*icebin_error)(-1, "run_timestep() not implemented");
    }

    /** Called at the beginning.  Returns variables in same place as run_timestep(),
    but doesn't actually run the timestep. */
    virtual void get_initial_state(double time_s)
    {
        (*icebin_error)(-1, "get_initial_state() not implemented");
    }
};      // class IceCoupler
// =========================================================
/** Serves as a base class for practical IceCouplers.  The
run_timestep() method is implemented, requiring subclasses to
implement the new method run_decoded().  Decoding converts a set of
(index, value) pairs into normal arrays (with NaN where no value was
given. */
class IceCoupler_Decode : public IceCoupler {
public :

    IceCoupler_Decode(IceCoupler::Type _type)
        : IceCoupler(_type) {}

    ~IceCoupler_Decode() {}

    /** @param index Index of each grid value.
    @param vals The values themselves -- could be SMB, Energy, something else...
    TODO: More params need to be added.  Time, return values, etc.
    @param time_s Time since start of simulation, in seconds */
    virtual void run_timestep(double time_s,
        blitz::Array<int,1> const &indices,
        std::vector<blitz::Array<double,1>> const &ivalsI);

    /** Runs a timestep after fields have been decoded.  This is what
    one will normally want to override, unless you wish to decode
    yourself. */
    virtual void run_decoded(double time_s,
        std::vector<blitz::Array<double,1>> const &ivalsI) = 0;
};

// =========================================================
#if 0
class IceCoupler_DISMAL : public IceCoupler_Decode
{

public:
    IceCoupler_DISMAL() : IceCoupler_Decode(IceCoupler::Type::DISMAL) {}

    /** @param index Index of each grid value.
    @param vals The values themselves -- could be SMB, Energy, something else...
    TODO: More params need to be added.  Time, return values, etc. */
    void run_decoded(double time_s,
        std::vector<blitz::Array<double,1>> const &valsI) {}
};
#endif
// =========================================================
class IceCoupler_Writer : public IceCoupler_Decode
{
    /** The IceCoupler we're affiliated with */
    IceCoupler const *main_model;

    /** Tells whether we are planning on writing the INPUT or OUTPUT fields */
    IceCoupler::IO io;

    // Dimensions to use when writing to netCDF
    std::vector<std::string> dim_names;
    std::vector<size_t> cur;        // Base index to write in netCDF
    std::vector<size_t> counts;

    // The output file we are writing to...
    std::string output_fname;

public:
    void init(IO _io, IceCoupler const *main_model);

    IceCoupler_Writer() :
        IceCoupler_Decode(IceCoupler::Type::WRITER) {}

    void start_time_set();

protected:
    bool output_file_initialized;
    /** This is called on-demand, the first time through run_decoded(). */
    void init_output_file();

public:
    /** @param index Index of each grid value.
    @param vals The values themselves -- could be SMB, Energy, something else...
    TODO: More params need to be added.  Time, return values, etc. */
    void run_decoded(double time_s,
        std::vector<blitz::Array<double,1>> const &ivalsI);

#if 0
protected:

    std::vector<netCDF::NcDim> add_dims(NcFile &nc);
    std::vector<netCDF::NcDim> get_dims(NcFile &nc);
#endif

};

}
