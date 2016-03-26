/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
class IceModel_Writer;

class IceModel {
public:
    BOOST_ENUM_VALUES( Type, int,
        (DISMAL)        (0)     // Demo Ice Sheet Model and LandIce
        (PISM)          (1)
        (ISSM)          (2)
        (WRITER)        (3)
    );
    const IceModel::Type type;

    /** Ordered specification of the variables (w/ units)
    to be passed IceBin->IceModel and IceModel->IceBin */
    enum IO {INPUT, OUTPUT};

//  friend class GCMCoupler;
//    friend class IceModel_Writer;
//  friend std::unique_ptr<IceModel> new_ice_model(IceModel::Type type,
//      GCMCoupler const *_coupler, IceRegridder const *_sheet);

public:
    GCMCoupler const *coupler;      // parent back-pointer
    IceRegridder *sheet;   // This is not const; see IceModel::update_ice_sheet()

    std::array<VarSet, 2> contract;     // [INPUT|OUTPUT]
    std::array<ibmisc::VarTransformer, 2> var_transformer;

    // Parameters provided by the GCM, to inform the coupling
    std::unique_ptr<GCMPerIceSheetParams> gcm_per_ice_sheet_params;

    // Writers called to record the input and output seen by this IceModel
    std::unique_ptr<IceModel_Writer> _iwriter, _owriter;

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

    /** @return the iwriter associated with this IceModel, if the
    IceModel is NOT responsible for calling the writer itself.
    If the IceModel WILL call the writer, returns NULL. */
    virtual IceModel_Writer *iwriter() { return _iwriter.get(); }
    virtual IceModel_Writer *owriter() { return _owriter.get(); }


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
    /** Allocate a new giss::CouplingContract, with the same lifetime as this IceModel. */
    VarSet *new_VarSet();
#endif

    IceModel(IceModel::Type _type) : type(_type) {}
    virtual ~IceModel();

    long ndata() const { return sheet->gridI->ndata(); }

    // --------------------------------------------------

public:

    /** Initialize any grid information, etc. from the IceSheet struct.
    @param vname_base Construct variable name from this, out of which to pull parameters from netCDF */
    virtual void ncread(ibmisc::NcIO &ncio, std::string const &vname_sheet);

    /** Event handler to let IceModels know the start time is (finally) set */
    virtual void start_time_set() {}

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

    /** Allows the IceModel to change the inputs used to create the
    regridding transformations.  This is used, for example, to make
    elev2 and mask2 consistent with an existing ice model (eg, PISM).
    It is called after init().
    Default implementation is to do nothing. */
    virtual void update_ice_sheet(ibmisc::NcIO &ncio, std::string const &vname) {}

};      // class IceModel
// =========================================================
/** Serves as a base class for practical IceModels.  The
run_timestep() method is implemented, requiring subclasses to
implement the new method run_decoded().  Decoding converts a set of
(index, value) pairs into normal arrays (with NaN where no value was
given. */
class IceModel_Decode : public IceModel {
public :

    IceModel_Decode(IceModel::Type _type)
        : IceModel(_type) {}

    ~IceModel_Decode() {}

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
class IceModel_DISMAL : public IceModel_Decode
{

public:
    IceModel_DISMAL() : IceModel_Decode(IceModel::Type::DISMAL) {}

    /** @param index Index of each grid value.
    @param vals The values themselves -- could be SMB, Energy, something else...
    TODO: More params need to be added.  Time, return values, etc. */
    void run_decoded(double time_s,
        std::vector<blitz::Array<double,1>> const &valsI) {}
};
#endif
// =========================================================
class IceModel_Writer : public IceModel_Decode
{
    /** The IceModel we're affiliated with */
    IceModel const *main_model;

    /** Tells whether we are planning on writing the INPUT or OUTPUT fields */
    IceModel::IO io;

    // Dimensions to use when writing to netCDF
    std::vector<std::string> dim_names;
    std::vector<size_t> cur;        // Base index to write in netCDF
    std::vector<size_t> counts;

    // The output file we are writing to...
    std::string output_fname;

public:
    void init(IO _io, IceModel const *main_model);

    IceModel_Writer() :
        IceModel_Decode(IceModel::Type::WRITER) {}

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


// ----------------------------------------------------------
struct SMBMsg {
    int sheetix;    // Ice sheet index, as defined by GCMRegridder::sheets_index
    int iI;         // Index of grid cell in ice sheet sheetix.
    double vals[1];     // Always at least one val; but this could be extended

    double &operator[](int i) { return *(vals + i); }

    /** @return size of the struct, given a certain number of values */
    static size_t size(int nfields)
        { return sizeof(SMBMsg) + (nfields-1) * sizeof(double); }

    static MPI_Datatype new_MPI_struct(int nfields);

    /** for use with qsort */
    static int compar(void const * a, void const * b);

};


class GCMCoupler {
public:
    /** Type tags for subclasses of GCMCoupler */
    BOOST_ENUM_VALUES( Type, int,
        (MODELE)        (0)
        (CESM)          (1)
    );
    Type const type;

    /** Filename this coupler (including grid) was read from. */
    std::string fname;
    /** Variable inside fname this coupler (including grid) was read from. */
    std::string vname;

    ibmisc::Domain<int> domainA;                // What's in our MPI halo?

    /** Main access to the core regridding of Icebin
    (for just this MPI node's part of the GCM domain) */
    GCMRegridder regridder;

    /** Access to entire regridding matrices, for all MPI nodes. */
    std::unique_ptr<GCMRegridder> regridder_full;

    /** Parameters passed from the GCM through to the ice model.
    These parameters cannot be specific to either the ice model or the GCM. */
    GCMParams gcm_params;

    /** See regridder.sheets_index */
    std::vector<std::unique_ptr<IceModel>> ice_models;


    ibmisc::UTSystem ut_system;     //!< Unit system for ConstantSets and CouplingContracts
    ibmisc::ConstantSet gcm_constants;      //!< Constants provided by the GCM

    /** Fields we receive from the GCM */
    VarSet gcm_outputs;

    /** Fields to send back to the GCM */
    VarSet gcm_inputs;


    /** Names of items used in the SCALARS dimension of VarTranslator.
    Used for ice_input and gcm_inputs.
    Scalars could be (eg) a timestep dt that is not known until runtime. */
    VarSet scalars;

    // Fields we read from the config file...

    /** File to which to write gcm_output.  If "", then don't write. */
    std::string gcm_out_file;

    /** File to which to write gcm_input.  (That is, stuff coming from
    Icebin back to the GCM.  If "", then don't write. */
    std::string gcm_in_file;

    GCMCoupler(Type _type) :
        type(_type),
        ut_system("")
    {
        gcm_constants.init(&ut_system);

        // Icebin requires orography on the ice grid, in order to
        // regrid in elevation space when things change.  Therefore, this
        // is added to the contract for all GCMs
        // gcm_inputs.add_field("elev2", "m", "ICE", "ice upper surface elevation");
        // No... this messes up what the GCM expects, and it's not used by the GCM.
        // Therefore, it should not be listed as a GCM input, it's a Icebin input.
        // Icebin will require it, somehow, as an IceModel output, and get it
        // directly from there.
    }

    virtual ~GCMCoupler() {}

    bool am_i_root() const { return (gcm_params.gcm_rank == gcm_params.gcm_root); }

    /** Read per-ice-sheet parameters that depend on the type of GCMCoupler. */
    virtual std::unique_ptr<GCMPerIceSheetParams>
    read_gcm_per_ice_sheet_params(
        ibmisc::NcIO &ncio,
        std::string const &sheet_vname) const = 0;

    virtual void ncread(
        std::string const &fname,
        std::string const &vname,
        ibmisc::Domain<int> &&domainA);

    void set_start_time(
        ibmisc::time::tm const &time_base,
        double time_start_s);


    /** Returns a unique rank number for each node in the parallel computation.
    Useful for debugging-type output. */
    int rank() const;

protected:
    /** @param time_s Time since start of simulation, in seconds
    Fields contained in the SMBMsg are INPUTS from the GCM.  They are
    therefore arranged according to gmc_inputs.  GCM inputs are converted
    into ice model inputs within IceModel::run_timestep(), which
    is called at the end of this method.
    @see gmc_inputs*/
    void call_ice_model(
        IceModel *model,
        int sheetno,
        double time_s,
        ibmisc::DynArray<SMBMsg> &rbuf,
        SMBMsg *begin, SMBMsg *end);


public:
    /** Top-level general-purpose method, called by icebin_modele.cpp
    (or other GCM-specific code).
    @param time_s Time (seconds) since the start of the GCM run.
    @param nfields Number of fields in sbuf.  Not all will necessarily be filled, in the case of heterogeneous ice models.
    @param sbuf the (filled) array of ice grid values for this MPI node.
    */
    void couple_to_ice(double time_s,
        int nfields,
        ibmisc::DynArray<SMBMsg> &sbuf,
        std::vector<SparseVector> &gcm_ivals);

    /** Follows the pattern of couple_to_ice()
    @param sbuf the (filled) array of ice grid values for this MPI node. */
    void get_initial_state(
    double time_s,
    std::vector<SparseVector> &gcm_ivals);  // Root node only: Already-allocated space to put output values.  Members as defined by the CouplingContract GCMCoupler::gcm_inputs

protected:
    void scaled_regrids(
        std::string const regrid_spec,
        std::vector<SparseMatrix> AvIs,
        SparseVector scaleA);

    void regrid_gcm_inputs_onroot(
        double time_s,
        std::vector<SparseVector> &gcm_ivals,   // Root node only: Already-allocated space to put output values.  Members as defined by the CouplingContract GCMCoupler::gcm_inputs
        unsigned int mask);

};

}
