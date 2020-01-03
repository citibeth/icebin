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
#include <icebin/e1ve0.hpp>

namespace icebin {

/** Different indices into gcm_inputs (see INDEXAE_* in LISheetIceBin.F90)*/
#if 0
BOOST_ENUM_VALUES( IndexAE, int,
    (A) (0)
    (E) (1)
    (ATOPO) (2)
    (ETOPO) (3)
)
#else
enum class IndexAE { A, E, ATOPO, ETOPO, COUNT};
static std::array<std::string, (unsigned)IndexAE::COUNT> const IndexAE_labels {"A", "E", "ATOPO", "ETOPO"};
static std::vector<char> indexae_grid {'A','E','A','E'};
#endif

template<int RANK>
    using TupleListLT = spsparse::TupleList<long,double,RANK>;

struct GCMInput {
    // http://www.boost.org/doc/libs/1_62_0/libs/serialization/doc/serialization.html#constructors
    friend class boost::serialization::access;

    // Contract inputs for the GCM on the A nad E grid, respectively (1D indexing).
    // _s = sparse indexing
    /** Contract inputs for the GCM on the A and E grid (See IndexAE) */
    std::vector<VectorMultivec> gcm_ivalss_s;

    // Regrid matrix to go from last step's elevation classes to this
    // step's elevation classes.
    // NOTE: Actual regrid matrix = I + E1vE0c
    spsparse::TupleList<int,double,2> E1vE0c;

    /** @param nvar Array specifying number of variables for each segment (A,E,ATOPO,ETOPO). */
    GCMInput(std::vector<int> const &nvar);
    /** @return Number of variables for each segment. */
    std::vector<int> nvar() const;

    void clear() {
        gcm_ivalss_s.clear();
        E1vE0c.clear();
    }

    template<class ArchiveT>
    void serialize(ArchiveT &ar, const unsigned int file_version)
    {
        ar & gcm_ivalss_s;
        ar & E1vE0c;
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

    int const icebin_base_hc = 0;    // First GCM elevation class that is an IceBin class (0-based indexing)

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

    // -------------- Info about the combined exchange grid (across all ice sheets)
    // This info is static.
    /** In combined exchange grid vector space: index of the start of
    each single-ice-sheet subspace. */
    std::vector<long> basesX;
    /** In combined exchange grid vector space: area of each cell */
    std::vector<double> areaX;


    // ------------ State that carries over across timesteps
    // Time period covered by the last time the coupler was called
    std::array<double,2> timespan; // {last_time_s, time_s}

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

    /** Description of fields to send back to the GCM; some VarSets are on E, some on A */
    std::vector<VarSet> gcm_inputs;
    /** Tells whether each gcm_inputs is on A or E grid */
    std::vector<char> gcm_inputs_grid;

    /** Names of items used in the SCALARS dimension of VarTranslator.
    Used for ice_input and gcm_inputs.
    Scalars could be (eg) a timestep dt that is not known until runtime. */
    VarSet scalars;

    /** XuE matrices from last timestep, used to compute E1vE0 */
    std::vector<std::unique_ptr<ibmisc::linear::Weighted_Eigen>> XuE0s;

    // Fields we read from the config file...

    GCMCoupler(Type _type, GCMParams &&_params);
    virtual ~GCMCoupler() {}

    /** Produces a date string in format YYMMDD */
    std::string sdate(double time_s) const;

    bool am_i_root() const { return gcm_params.am_i_root(); }

    /** Run the coupling procedure for all ice sheets.
    @return All output from copuling to pass back to the GCM */
    virtual GCMInput couple(
        double time_s,        // Simulation time [s]
        VectorMultivec const &gcm_ovalsE,
        bool run_ice);    // if false, only initialize

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
