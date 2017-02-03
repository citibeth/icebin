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

#include <ibmisc/VarTransformer.hpp>
#include <ibmisc/DynArray.hpp>
#include <ibmisc/udunits2.hpp>
#include <ibmisc/ConstantSet.hpp>

#include <icebin/IceCoupler.hpp>
#include <icebin/GCMParams.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/VarSet.hpp>
#include <icebin/multivec.hpp>

namespace icebin {


template<int RANK>
    using TupleListLT = spsparse::TupleList<long,double,RANK>;

struct GCMInput {
    // http://www.boost.org/doc/libs/1_62_0/libs/serialization/doc/serialization.html#constructors
    friend class boost::serialization::access;

    // Mapping from the index of a variable in gcm_ivalsE/gcm_ivalsA
    // and the index within the GCMCoupler::gcm_inputs
    std::array<VectorMultivec, GridAE::count> gcm_ivalsAE;

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
    TupleListLT<2> E1vE0_s;

    // Regrid matrix to convert to atmosphere.
    // (EvA is assumed by GCM, as long as AvE is local; see Fischer&Nowicki 2014)
    TupleListLT<2> AvE1_s;
    TupleListLT<1> wAvE1_s;

    // Used for temperature downscaling according to a lapse rate
    TupleListLT<1> elevE1_s;

    GCMInput(std::array<int, GridAE::count> const &nvar) :
        gcm_ivalsAE({
            VectorMultivec(nvar[0]),
            VectorMultivec(nvar[1]),
        })
    {}

    std::array<int, GridAE::count> nvar() const
    {
        std::array<int, GridAE::count> ret;
        for (int i=0; i<GridAE::count; ++i) ret[i] = gcm_ivalsAE[i].nvar;
        return ret;
    }

    template<class ArchiveT>
    void serialize(ArchiveT &ar, const unsigned int file_version)
    {
        ar & gcm_ivalsAE;
        ar & E1vE0_s;
        ar & AvE1_s;
        ar & wAvE1_s;
        ar & elevE1_s;
    }
};

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

    /** Main access to the core regridding of Icebin */
    GCMRegridder gcm_regridder;

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
protected:
    virtual int _read_nhc_gcm() = 0;

public:
    /** See regridder.sheets_index */
    std::vector<std::unique_ptr<IceCoupler>> ice_couplers;

    ibmisc::UTSystem ut_system;     //!< Unit system for ConstantSets and CouplingContracts
    ibmisc::ConstantSet gcm_constants;      //!< Constants provided by the GCM

    /** Description of fields we receive from the GCM; all on the E grid. */
    VarSet gcm_outputsE;

    /** Description of fields to send back to the GCM; some on E, some on A */
    std::array<VarSet, GridAE::count> gcm_inputsAE;    // gcm_inputsE, gcm_inputsA

    /** Names of items used in the SCALARS dimension of VarTranslator.
    Used for ice_input and gcm_inputs.
    Scalars could be (eg) a timestep dt that is not known until runtime. */
    VarSet scalars;

    // Fields we read from the config file...

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
        // Icebin will require it, somehow, as an IceCoupler output, and get it
        // directly from there.
    }

    virtual ~GCMCoupler() {}

    virtual void ncread(
        std::string const &grid_fname,
        std::string const &config_fname,
        std::string const &vname);

    /** Private; called from gcmce_cold_start() */
    void cold_start(
        ibmisc::Datetime _time_base,
        double time_start_s);

    /** @param am_i_root
        Call with true if calling from MPI root; false otherwise.
        The core coupling/regridding computation only runs on root.
        But other MPI ranks need to go along for the ride, assuming that
        the ice model uses MPI. */
    GCMInput couple(
        double time_s,        // Simulation time [s]
        VectorMultivec const &gcm_ovalsE,
        bool run_ice,
        bool am_i_root);

    /** Top-level ncio() to log output from coupler. (coupler->GCM) */
    void ncio_gcm_input(ibmisc::NcIO &ncio,
        GCMInput &out,        // All MPI ranks included here
        std::array<double,2> &timespan,    // timespan[1] = current time_s
        std::string const &time_units,
        std::string const &vname_base);

    /** Top-level ncio() to log input to coupler. (GCM->coupler) */
    void ncio_gcm_output(ibmisc::NcIO &ncio,
        VectorMultivec const &gcm_ovalsE,
        std::array<double,2> &timespan,    // timespan[1] = current time_s
        std::string const &time_units,
        std::string const &vname_base);

};

}
