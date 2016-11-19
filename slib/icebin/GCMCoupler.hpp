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
#include <icebin/GCMPerIceSheetParams.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/VarSet.hpp>

namespace icebin {

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
    std::vector<std::unique_ptr<IceCoupler>> ice_couplers;


    ibmisc::UTSystem ut_system;     //!< Unit system for ConstantSets and CouplingContracts
    ibmisc::ConstantSet gcm_constants;      //!< Constants provided by the GCM

    /** Fields we receive from the GCM; all on the E grid. */
    VarSet gcm_outputsE;

    /** Fields to send back to the GCM; some on E, some on A */
    enum class GCMI { E, A, COUNT };
    VarSet[GCMI::COUNT] gcm_inputs;    // gcm_inputsE, gcm_inputsA


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
        // Icebin will require it, somehow, as an IceCoupler output, and get it
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
    into ice model inputs within IceCoupler::run_timestep(), which
    is called at the end of this method.
    @see gmc_inputs*/
    void call_ice_model(
        IceCoupler *model,
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
