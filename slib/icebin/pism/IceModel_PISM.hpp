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

// --------------------------------
// PISM Includes... want to be included first
#include <petsc.h>
#include <petscvec.h>
#include <pism/base/util/Context.hh>
#include <pism/base/util/IceGrid.hh>
#include <pism/base/util/iceModelVec.hh>
#include <pism/base/iceModel.hh>

#include <pism/base/util/pism_options.hh>
#include <pism/coupler/atmosphere/PAFactory.hh>
#include <pism/coupler/ocean/POFactory.hh>
#include <pism/coupler/surface/PSFactory.hh>

#include <pism/base/util/PISMTime.hh>

#include <pism/icebin/IBIceModel.hh>
#include <pism/icebin/VecBundleWriter.hh>
#include <pism/icebin/IBSurfaceModel.hh>
#include <pism/base/util/petscwrappers/PetscInitializer.hh>
// --------------------------------
#include <mpi.h>
#include <icebin/GCMCoupler.hpp>
#include <memory>
#include <icebin/Grid.hpp>
#include <icebin/GCMParams.hpp>


namespace icebin {
namespace gpism {

// =================================================


class IceModel_PISM : public IceModel
{
    icebin::Grid_XY const *icebin_gridI;
    MPI_Comm pism_comm;         // Commnicator used by ice model
    PetscMPIInt _pism_rank, _pism_size;
    const int pism_root = 0;

public:
    PetscMPIInt pism_rank() const { return _pism_rank; }
    PetscMPIInt pism_size() const { return _pism_size; }
    bool am_i_root() const { return _pism_rank == 0; }

private:
    /** Arguments to be passed to PISM at PISM initialization */
    std::vector<std::string> pism_args;

    /** Should we write inputs to PISM via PISM's dump() functionality?
    This is (almost certainly) supercede by the IceModel_Writer class,
    except in cases of extreme debugging. */
    const bool write_pism_inputs = false;

    /** Constants used by PISM, and fetched from GCM config */

    // These should be declared in the same order they're created,
    // so they get destroyed in the proper reverse order.
    // See: http://msdn.microsoft.com/en-us/library/8183zf3x%28v=vs.110%29.aspx
    std::unique_ptr<pism::petsc::Initializer> petsc_initializer;
    pism::IceGrid::Ptr pism_grid;
    std::unique_ptr<pism::icebin::IBIceModel> pism_ice_model;
    pism::icebin::IBSurfaceModel *pism_surface_model;   // We don't own this.
public:
    pism::Context::ConfigPtr pism_config() { return pism_ice_model->ctx()->config(); }

private:
    // Stuff used for Scatter/Gather
    pism::petsc::DM::Ptr da2;
    Vec g2, g2natural;  //!< global Vecs used to transfer data to/from processor 0.
    VecScatter scatter; //!< VecScatter used to transfer data to/from processor 0.
    pism::petsc::Vec::Ptr Hp0;            //!< Resulting vector on process 0

    // Corresponding PISM variable for each input field
    std::vector<pism::IceModelVec2S *> pism_ivars;

    // --- Corresponding PISM variable for each output field
    // The variable, as it exists in PISM
    std::vector<pism::IceModelVec2S const *> pism_ovars;

//  // An MPI-collected version of the variable
//  std::vector<blitz::Array<double,2> icebin_ovars;


    double BY_ICE_DENSITY;      // CONSTANT Used to prepare input for PISM

    /** Should we upate the elevation field in update_ice_sheet()?  Normally, yes.
    But in some TEST CASES ONLY --- when the SMB field was created with a different
    set of elevations than the ice model is using --- then this can cause problems
    in the generated SMB fields. */
    bool update_elevation = true;

    // NetCDF output files
    std::unique_ptr<pism::icebin::VecBundleWriter> pism_in_nc, pism_out_nc;

    // ------------------------
public:

    virtual void ncread(ibmisc::NcIO &ncio, std::string const &vname_sheet);

    int nx() { return pism_grid->Mx(); }
    int ny() { return pism_grid->My(); }

    // ----------------------------------------------

    // ----------------------------------------------
    void iceModelVec2S_to_blitz_xy(pism::IceModelVec2S const &pism_var,
        blitz::Array<double,2> &ret);
    // ----------------------------------------------

    void update_ice_sheet(ibmisc::NcIO &ncio, std::string const &vname);

    IceModel_PISM();

    virtual ~IceModel_PISM()
        { deallocate(); }

protected:

    void allocate();

    void deallocate();

public:
    /** Event handler to let IceModels know the start time is (finally) set */
    void start_time_set()
        { allocate(); }

    void setup_contracts_modele();

    /** Transfers a constant from GCMCoupler::gcm_constants to PISM's configuration variable.
    Runs from within transfer_constants_xxxx() */
    void transfer_constant(std::string const &dest, std::string const &src, double multiply_by=1.0, bool set_new = false);

    /** @param set_new If true, PISM constant will be set, even if it
    was not already set in the configuration.  This defaults to false,
    as a program error check against misspelled parameter names. */
    void set_constant(std::string const &dest, double src_val, std::string const &src_units, bool set_new = false);

public:
    void run_timestep(double time_s,
        blitz::Array<int,1> const &indices,
        std::vector<blitz::Array<double,1>> const &ivals2);

    void get_initial_state(double time_s);

private:

    /** Copies PISM->Icebin output variables from PISM variables to
    the Icebin-supplied variables (on the root node).
    @param mask Only do it for variables where (flags & mask) == mask.  Set to 0 for "all." */
    void get_state(unsigned int mask = 0);
};

}}  // namespace icebin::pism
