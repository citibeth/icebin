/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <cstdlib>

#include <boost/mpi.hpp>

#include <ibmisc/VarTransformer.hpp>
#include <ibmisc/DynArray.hpp>
#include <ibmisc/udunits2.hpp>
#include <ibmisc/ConstantSet.hpp>

#include <icebin/IceCoupler.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/VarSet.hpp>
#include <icebin/multivec.hpp>

namespace icebin {


template<int RANK>
    using TupleListLT = spsparse::TupleList<long,double,RANK>;

struct GCMInput {
    // http://www.boost.org/doc/libs/1_62_0/libs/serialization/doc/serialization.html#constructors
    friend class boost::serialization::access;

    // Contract inputs for the GCM on the A nad E grid, respectively (1D indexing).
    // _s = sparse indexing
    std::vector<VectorMultivec> gcm_ivalss_s;
    std::vector<std::vector<double>> gcm_ivalss_weight_s;

    // Output fields from coupling back into GCM
    //    * Computed by merging from IceCoupler::couple()
    //    * UNSCALED: values must be divided by last column (the weight)
//    VectorMultivec gcm_ivalsA_s;
//    VectorMultivec gcm_ivalsE_s;


    // This is all for elevation classes in ICE space (ice_nhc, not gcm_nhc)

    // Values required to update TOPO, etc. in ModelE
    // (see add_fhc.py for how these are to be used)
    // We can get these from AvE
    // TupleList<1> wAvE;     // Area of A grid cells that overlap ice
    //TupleList<1> areaA;    // Total (native) area of A grid cells
    //TupleList<1> elevA;    // Used for atmosphere orography (ZATMO in ModelE)


    // _s means these matrices and vectors all use sparse indexing

    // Regrid matrix to go from last step's elevation classes to this
    // step's elevation classes.
    linear::Weighted_Tuple E1vE0_unscaled;    // Sparse indexing

    // Regrid matrix to convert to atmosphere.
    // (EvA is assumed by GCM, as long as AvE is local; see Fischer&Nowicki 2014)
//    linear::Weighted_Tuple AvE1_unscaled;

Root arrays should be returned...


    GCMInput(std::vector<int> const &nvar) {
        for (int nv : nvar) gcm_ivals_s.push_back(VectorMultivec(nv));
    }

    std::array<int, GridAE::count> nvar() const
    {
        std::array<int, GridAE::count> ret;
        for (int i=0; i<GridAE::count; ++i) ret[i] = gcm_ivalsAE_s[i].nvar;
        return ret;
    }

    template<class ArchiveT>
    void serialize(ArchiveT &ar, const unsigned int file_version)
    {
        ar & gcm_ivalsAE_s;
        ar & E1vE0_s;
        ar & AvE1_s;
        ar & wAvE1_s;
    }
};
// =============================================================================
/** A segment of elevation classes (see add_fhc.py) */
struct HCSegmentData {
    std::string name;
    int base;    // First elevation class of this segment
    int size;    // Number of elevation classes in this segment

    HCSegmentData(std::string const &_name, int _base, int _size)
        : name(_name), base(_base), size(_size) {}
};

extern HCSegmentData &get_segment(std::vector<HCSegmentData> &hc_segments, std::string const &name);

inline HCSegmentData const &get_segment(std::vector<HCSegmentData> const &hc_segments, std::string const &name)
    { return get_segment(const_cast<std::vector<HCSegmentData> &>(hc_segments), name); }

/** Parses a spec. string (eg: "legacy,sealand,ec") to a usable set of HCSegments. */
extern std::vector<HCSegmentData> parse_hc_segments(std::string const &str);

/** Parameters passed from the GCM through to the ice model.
These parameters cannot be specific to either the ice model or the GCM.
TODO: Make procedure to read rundeck params and set this stuff up. */
struct GCMParams {
    // ------- Initialized in constructor
    MPI_Comm gcm_comm;
    int gcm_root;        // Root of the MPI group
    int gcm_rank;            // MPI rank of this node
    bool am_i_root() const { return gcm_rank == gcm_root; }
    boost::mpi::communicator world;

    // ------- Passed into GCMCoupler::allocate()

    // Name of the IceBin config file.
    // Other input files should be listed here as abosolute paths.
    std::string icebin_config_fname;

    bool icebin_logging = true ;    // Should IceBin log input & output?

    // Should IceBin update topography?
    bool dynamic_topo = false;

    std::vector<HCSegmentData> hc_segments {    // 0-based
        HCSegmentData("legacy", 0, 1),
        HCSegmentData("sealand", 1, 2),
        HCSegmentData("ec", 3, -1)};    // Last segment must be called ec
    int icebin_base_hc;    // First GCM elevation class that is an IceBin class (0-based indexing)
    std::string primary_segment = "ec";

    HCSegmentData &segment(std::string const &name)
        { return get_segment(hc_segments, name); }
    HCSegmentData const &segment(std::string const &name) const
        { return get_segment(hc_segments, name); }

    GCMParams(MPI_Comm _gcm_comm, int _gcm_root);
};
// =============================================================================
class GCMCoupler {
public:
    /** Type tags for subclasses of GCMCoupler */
    BOOST_ENUM_VALUES( Type, int,
        (MODELE)        (0)
        (CESM)          (1)
    );
    Type const type;

    // ------- Set in GCMCoupler::cold_start()
    ibmisc::Datetime time_base;    // yy,mm,dd,hh,mm,ss
    ibmisc::TimeUnit time_unit;    // Equiv. to CF-compliant time unit string
    double time_start_s;        // Start of simulation, as far as ice model is concerned (seconds since time_base).


    // Last time this coupler was called
    double last_time_s;

    /** Filename this coupler (including grid) was read from. */
    std::string icebin_in;

    /** Parameters read from IceBin config file */
    std::string output_dir;

    /** Set to false and IceBin will pass a zero SMB and appropriately
    zero B.C. to the ice sheet.  This is for testing. */
    bool use_smb;

    /** Main access to the core regridding of Icebin */
    std::shared_ptr<GCMRegridder> gcm_regridder;

    /** Parameters (not physical constants) passed from the GCM
    through to the ice model.  These parameters cannot be specific to
    either the ice model or the GCM. */
    GCMParams gcm_params;

    /** Number of elevation classes the GCM sees */
    int _nhc_gcm = -1;
    int nhc_gcm() {
        if (_nhc_gcm < 0) _nhc_gcm = _read_nhc_gcm();
        return _nhc_gcm;
    }

    long nE_gcm()
        { return gcm_regridder->nA() * nhc_gcm(); }

protected:
    virtual int _read_nhc_gcm() = 0;

public:

    /** See regridder.ice_regridders().index */
    std::vector<std::unique_ptr<IceCoupler>> ice_couplers;

    ibmisc::UTSystem ut_system;     //!< Unit system for ConstantSets and CouplingContracts
    ibmisc::ConstantSet gcm_constants;      //!< Constants provided by the GCM

    /** Description of fields we receive from the GCM; all on the E grid. */
    VarSet gcm_outputsE;

    /** Description of fields to send back to the GCM; some on E, some on A */
    std::vector<VarSet> gcm_inputs;
    /** Tells whether each gcm_inputs is on A or E grid */
    std::vector<char> gcm_inputs_grid;

    /** Names of items used in the SCALARS dimension of VarTranslator.
    Used for ice_input and gcm_inputs.
    Scalars could be (eg) a timestep dt that is not known until runtime. */
    VarSet scalars;

    // Fields we read from the config file...

    GCMCoupler(Type _type, GCMParams &&_params);
    virtual ~GCMCoupler() {}

    /** Produces a date string in format YYMMDD */
    std::string sdate(double time_s) const;

    bool am_i_root() const { return gcm_params.am_i_root(); }

protected:
    virtual void _ncread(
        ibmisc::NcIO &ncio_config,
        std::string const &vname);

public:
    void ncread(
        std::string const &config_fname,
        std::string const &vname)
    {
        ibmisc::NcIO ncio(config_fname, netCDF::NcFile::read);
        _ncread(ncio, vname);
    }

    /** Locates an ice model input file, according to resolution rules
        of the GCM. */
    virtual std::string locate_input_file(
        std::string const &sheet_name,        // eg: greenland
        std::string const &file_name) = 0;        // eg: pism_Greenland_5km_v1.1.nc

    /** Private; called from gcmce_cold_start() */
    void cold_start(
        ibmisc::Datetime _time_base,
        double time_start_s);

    /** Top level method to re-compute values originally loaded from the TOPO file. */
    virtual void update_topo(double time_s) = 0;

    /** @param am_i_root
        Call with true if calling from MPI root; false otherwise.
        The core coupling/regridding computation only runs on root.
        But other MPI ranks need to go along for the ride, assuming that
        the ice model uses MPI. */
    GCMInput couple(
        double time_s,        // Simulation time [s]
        VectorMultivec const &gcm_ovalsE,
        bool run_ice);

    /** Top-level ncio() to log output from coupler. (coupler->GCM) */
    void ncio_gcm_input(ibmisc::NcIO &ncio,
        GCMInput &out,        // All MPI ranks included here
        std::array<double,2> &timespan,    // timespan[1] = current time_s
        ibmisc::TimeUnit const &time_unit,
        std::string const &vname_base);

    /** Top-level ncio() to log input to coupler. (GCM->coupler) */
    void ncio_gcm_output(ibmisc::NcIO &ncio,
        VectorMultivec const &gcm_ovalsE,
        std::array<double,2> &timespan,    // timespan[1] = current time_s
        ibmisc::TimeUnit const &time_unit,
        std::string const &vname_base);

};

}
