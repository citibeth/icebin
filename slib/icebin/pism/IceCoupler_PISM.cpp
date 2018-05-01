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

#include <mpi.h>    // For Intel MPI, mpi.h must be included before stdio.h

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <pism/base/stressbalance/PISMStressBalance.hh>
#include <pism/earth/PISMBedDef.hh>

#include <ibmisc/netcdf.hpp>
#include <ibmisc/ibmisc.hpp>

#include <icebin/GCMCoupler.hpp>
#include <icebin/pism/IceCoupler_PISM.hpp>
#include <icebin/contracts/contracts.hpp>

extern "C" void libpismutil_refaddr();

using namespace ibmisc;
using namespace pism;
using namespace netCDF;

namespace icebin {
namespace gpism {

static double const nan = std::numeric_limits<double>::quiet_NaN();

IceCoupler_PISM::IceCoupler_PISM()
    : IceCoupler(IceCoupler::Type::PISM),
    write_pism_inputs(true)
{
}

#if 0
// PISM command line arguments that are paths, and thus need pathname
// resolution.
// For PISM stable0.5 branch
// static std::set<std::string> path_args = {"config_override", "i", "o", "surface_given_file", "extra_file", "ts_file"};
// For dev branch
static std::map<std::string, std::string> path_args = {
    {"i", "i"},
    {"surface_given_file", "i"},
    {"ocean_kill", "i"},
    {"ocean_kill_file", "i"},
    {"o", "o"},
    {"extra_file", "o"},
    {"ts_file", "o"}};
#endif

void IceCoupler_PISM::ncread(ibmisc::NcIO &ncio_config, std::string const &vname_sheet)
{
    IceCoupler::ncread(ncio_config, vname_sheet);

    printf("BEGIN IceCoupler_PISM::ncread(%s)\n", vname_sheet.c_str());
    GCMParams const &_gcm_params(gcm_coupler->gcm_params);

    this->icebin_specI = ice_regridder ?
        dynamic_cast<GridSpec_XY const *>(&*agridI().spec)
        : NULL;

    // General args passed to the ice sheet, regardless of which ice model is being used
    NcVar info_var(ncio_config.nc->getVar(vname_sheet + ".info"));
    // PISM parameters, passed to PISM via argv
    NcVar pism_var(ncio_config.nc->getVar(vname_sheet + ".pism"));

    // Get simple arguments
    get_or_put_att(info_var, 'r', "update_elevation", &update_elevation, 1);
    get_or_put_att(info_var, 'r', "output_dir", output_dir);

    // Create arguments from PISM configuration
    pism_args.push_back("icebin_pism");

    // Get arguments from IceBin configuration
    std::map<std::string, NcVarAtt> pism_atts(pism_var.getAtts());
    for (auto jj=pism_atts.begin(); jj != pism_atts.end(); ++jj) {
        std::string const &name(jj->first);
        NcAtt const &att(jj->second);
        std::string val;
        att.getValues(val);

        pism_args.push_back("-" + name);
        pism_args.push_back(val);
    }

#if 0
    // Hard-code these variables because we will always need them
    pism_args.push_back("-extra_vars");
    pism_args.push_back("climatic_mass_balance_cumulative,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,flux_divergence");
#endif
    printf("END IceCoupler_PISM::ncread()\n");
}
// ======================================================================
void IceCoupler_PISM::_cold_start(
    ibmisc::Datetime const &time_base,
    double time_start_s)
{
    printf("BEGIN IceCouple_PISM::_cold_start()\n");

    // ------- Now instantiate PISM!
    // Convert PISM arguments to old C style
    int argc = pism_args.size();
    char *argv_array[argc];
    std::vector<char> all_str;
    for (int i=0; i<argc; ++i) {
        std::string &arg = pism_args[i];
        for (unsigned int j=0; j<arg.size(); ++j) all_str.push_back(arg[j]);
        all_str.push_back('\0');
    }
    char *pos = &all_str[0];
    for (int i=0; i<argc; ++i) {
        std::string &arg = pism_args[i];
        argv_array[i] = pos;
        pos +=arg.size() + 1;
    }
    char **argv = argv_array;

printf("*** PISM Args:");
for (int i=0; i<argc; ++i) {
    printf(" %s", argv[i]);
}
printf("\n");


    // Set up communicator for PISM to use
    // Use same group of processes.
    // No spawning or intercommunicators for now --- maybe not ever.
//  MPI_Comm_dup(gcm_params.gcm_comm, &pism_comm);
    pism_comm = gcm_coupler->gcm_params.gcm_comm;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(pism_comm, &_pism_rank); PISM_CHK(ierr, "MPI_Comm_rank");
    ierr = MPI_Comm_size(pism_comm, &_pism_size); PISM_CHK(ierr, "MPI_Comm_size");

printf("[%d] pism_size = %d\n", pism_rank(), pism_size());

    // -------------- Initialize Petsc
    // From personal correspondence with Barry Smith <bsmith@mcs.anl.gov>
    //  > 3. The way PETSc provides stacktraces is clever.  However, it
    //   doesn't really work well with the non-PETSc parts of the program
    //   which don't (and won't) follow PETSc conventions.  It's also
    //   clumsy to use, i.e. checking error returns on
    //   every. single. call.  And finally, there are easier and less
    //   invasive ways to get stacktraces (which I am discovering and
    //   implementing).  Again... stacktraces are nice, but they have
    //   nothing to do with the core functionality of PETSc in this case.
    //
    // You can have PETSc either abort on errors with either
    // PetscOptionsSetValue("-on_error_abort","true") or
    // PetscOptionsSetValue("-on_error_mpiabort","true") BEFORE
    // PetscInitialize() or you can call PetscPushErrorHandler() after
    // PetscInitialize() to set a custom or standard error handler, for
    // example PetscReturnErrorHandler() simply returns the error flag up
    // the stack without printing the stack
    // trace. PetscAbortErrorHandler() and PetscMPIAbortErrorHandler()
    // trigger an immediate abort or MPI_Abort().  If you don't want to
    // check error codes always in C++ you could possible generate an
    // exception on the errors and optionally catch them on the way up
    // but I don't know how to a handle that in Fortran, plus at the
    // interface calls between C++ and C or Fortran you would need to
    // catch all the exceptions and translate them to error codes again.
    //
    // Also see here for proper calling sequence of PetscOptionsSetValue()
    //    http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscOptionsSetValue.html
    // petsc_initializer.reset(new pism::petsc::Initializer(pism_comm, argc, argv));
//    printf("Setting -on_error_mpiabort\n");
//    PetscOptionsSetValue(NULL, "-on_error_mpiabort", "true");
    printf("Doing -no_signal_handler\n");
    PetscOptionsSetValue(NULL, "-no_signal_handler", "true");
    petsc_initializer.reset(new pism::petsc::Initializer(argc, argv, "IceBin GCM Coupler"));
//    printf("Calling PetscPopErrorHandler()\n");
//    PetscPopErrorHandler();
    // ------------------------------------

    // verbosityLevelFromOptions();    // https://github.com/pism/pism/commit/3c75fd63
    Context::Ptr ctx = context_from_options(pism_comm, "IceCoupler_PISM");
    Logger::Ptr log = ctx->log();

    log->message(2, "IceBin %s (GCM Coupler)\n",
                 PISM_Revision);

    bool input_file_set = options::Bool("-i", "input file name");
    std::string usage =
      "  pismr -i IN.nc [-bootstrap] [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i          IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -bootstrap  enable heuristics to produce an initial state from an incomplete input\n"
      "notes:\n"
      "  * option -i is required\n"
      "  * if -bootstrap is used then also '-Mx A -My B -Mz C -Lz D' are required\n";


    if (not input_file_set) {
      (*icebin_error)(-1, "PISM option -i is required");
    } else {
      std::vector<std::string> required;
      required.clear();

      if (show_usage_check_req_opts(*log, "pismr", required, usage)) {
          (*icebin_error)(-1, "Invalid PISM options");
      }
    }

    options::String profiling_log = options::String("-profile",
        "Save detailed profiling data to a file.");

    Config::Ptr config = ctx->config();

    // ------------------------------ \\
    // Get arguments from the GCM
    ibmisc::Datetime const &tb(gcm_coupler->time_base);
    std::string reference_date = (boost::format("%04d-%02d-%02d") % tb.year() % tb.month() % tb.day()).str();
    config->set_string("time.reference_date", reference_date);
    // ------------------------------ //

#if 0    // Don't bother with profiling inside of GCM
    if (profiling_log.is_set()) {
      ctx->profiling().start();
    }
#endif

    log->message(3, "* Setting the computational grid...\n");
    pism_grid = IceGrid::FromOptions(ctx);

    pism::icebin::IBIceModel::Params params;
        params.time_start_s = gcm_coupler->time_start_s;

        params.output_dir = output_dir;

    pism_ice_model.reset(new pism::icebin::IBIceModel(pism_grid, ctx, params));

    if (nx() * ny() != ice_regridder->nI()) (*icebin_error)(-1,
        "nI does not match: %ld vs. %ld\n", nx()*ny(), ice_regridder->nI());

    // ------------------------------------------- \\

    // Transfer constants from GCM to PISM, and also set up coupling contracts.
    // This is the right place to do it, since PISM is now fully up and functional,
    // and all PISM config files have been read.
    // This call through the GCMCoupler will call back to setup_contracts_xxx().
    contracts::setup(gcm_coupler, this);

//printf("start = %f\n", pism_grid->time->start());
//printf("end = %f\n", pism_grid->time->end());

    // This has the following stack trace:
    //  IceCoupler::init()                    [iceModel.cc]
    //  IceCoupler::model_state_setup()       [iMinit.cc]
    //  IceCoupler::init_couplers()           [iMinit.cc]
    //  surface->init()
    pism_ice_model->init();

    // ============== Set up variables for INPUT contract

    // During the pism_ice_model->init() call above the PISMIceModel
    // class (derived from PISM's IceModel) allocated an instance
    // of PSConstantICEBIN. This instance is owned and will be
    // de-allocated by PISMIceModel pism_ice_model.
    pism_surface_model = pism_ice_model->ib_surface_model();

    // Set up corresponence between IceBin fields and variables
    // in the PISM data structures.
    int ix;
    pism_ivars.resize(contract[INPUT].size(), NULL);

    // We don't really use this, but we do need to store and pass through for conservation computations
    ix = contract[INPUT].index.at("massxfer");        // [kg m-2 s-1]
        pism_ivars[ix] = &pism_surface_model->massxfer;

    // Ignore surface_temp, it is not useful...
    ix = contract[INPUT].index.at("enthxfer");
        pism_ivars[ix] = &pism_surface_model->enthxfer;

    ix = contract[INPUT].index.at("ice_top_bc_temp");
        pism_ivars[ix] = &pism_surface_model->ice_top_bc_temp;
    ix = contract[INPUT].index.at("ice_top_bc_wc");
        pism_ivars[ix] = &pism_surface_model->ice_top_bc_wc;



    // Check that all PISM inputs are bound to a variable
    bool err = false;
    for (unsigned i=0; i<pism_ivars.size(); ++i) {
        IceModelVec2S *pism_var = pism_ivars[i];
        if (!pism_var && !(contract[INPUT][i].flags & contracts::PRIVATE)) {
            fprintf(stderr,
                "PISM input %s (i=%d) is not bound to a variable\n",
                contract[INPUT].data[i].name.c_str(), i);
            err = true;
        }
    }
    if (err) (*icebin_error)(-1, "Exiting due to errors");


    // Initialize scatter/gather stuff
    vtmp.create(pism_grid, "vtmp", WITHOUT_GHOSTS);
    vtmp_p0 = vtmp.allocate_proc0_copy();

    // Initialize scatter/gather stuff (probably obsolete)
    da2 = pism_grid->get_dm(1, // dm_dof
        pism_grid->ctx()->config()->get_double("grid.max_stencil_width"));

    ierr = DMCreateGlobalVector(*da2, &g2); PISM_CHK(ierr, "DMCreateGlobalVector");

    // note we want a global Vec but reordered in the natural ordering
    // so when it is scattered to proc zero it is not all messed up;
    // see above
    ierr = DMDACreateNaturalVector(*da2, &g2natural);
        PISM_CHK(ierr, "DMDACreateNaturalVector");

//    // next get context *and* allocate samplep0 (on proc zero only, naturally)
//    ierr = VecScatterCreateToZero(g2natural, &scatter, &Hp0);
//        PISM_CHK(ierr, "VecScatterCreateToZero");


    // ============== Set up variables for OUTPUT contract
    // -------------- Allocate blitz:Array<double,1> output variables,
    // which are passed back to Icebin.

    // -------------- Link to PISM-format output variables, used to fill ovars
    pism_ovars.resize(contract[OUTPUT].size(), NULL);
    ix = contract[OUTPUT].index.at("ice_top_elevation");       // Elevation of top surface of ice sheet
        pism_ovars[ix] = &pism_ice_model->ice_surface_elevation(); // see PISM's iceModel.hh
    ix = contract[OUTPUT].index.at("ice_thickness");
        pism_ovars[ix] = &pism_ice_model->ice_thickness();
    ix = contract[OUTPUT].index.at("bed_topography");
        pism_ovars[ix] = &pism_ice_model->bed_model()->bed_elevation();

    ix = contract[OUTPUT].index.at("mask");
        pism_ovars[ix] = &pism_ice_model->cell_type();
    ix = contract[OUTPUT].index.at("elev");
        pism_ovars[ix] = &pism_ice_model->elevmaskI;



    // Specific enthalpy (PISM-style) of the top layer
    ix = contract[OUTPUT].index.at("ice_top_senth");
        pism_ovars[ix] = &pism_ice_model->ice_top_senth;

    // For MassEnergyBudget variables that have a contract name specified,
    // link them up into pism_ovars now.
    for (auto ii = pism_ice_model->rate.all_vecs.begin(); ii != pism_ice_model->rate.all_vecs.end(); ++ii) {
        if (ii->contract_name == "") continue;

        int ix = contract[OUTPUT].index.at(ii->contract_name);
        pism_ovars[ix] = &ii->vec;
    }

    // -------------- Initialize pism-out.nc
    {
        boost::filesystem::path output_dir(params.output_dir);
        std::string ofname = (output_dir / "pism-out.nc").string();

        // Convert array from pism::IceModelVec2S* to pism::IceModelVec*
        std::vector<pism::IceModelVec const *> vecs;
        for (unsigned i=0; i<pism_ovars.size(); ++i) {
            if (!(contract[OUTPUT][i].flags & contracts::PRIVATE))
                vecs.push_back(pism_ovars[i]);
        }

        pism_out_nc.reset(new pism::icebin::VecBundleWriter(
            pism_ice_model->grid(), ofname, vecs));
        pism_out_nc->init();
    }

    // ------------- Initialize pism-in.nc
    {
        boost::filesystem::path output_dir(params.output_dir);
        std::string ofname = (output_dir / "pism-in.nc").string();
        std::vector<pism::IceModelVec const *> vecs;
        for (unsigned i=0; i<pism_ivars.size(); ++i) {
            if (!(contract[INPUT][i].flags & contracts::PRIVATE))
                vecs.push_back(pism_ivars[i]);
        }
        pism_in_nc.reset(new pism::icebin::VecBundleWriter(
            pism_ice_model->grid(), ofname, vecs));
        pism_in_nc->init();
    }


    // ============== Miscellaneous
    // Check that grid dimensions match
    if (icebin_specI) {
        if ((pism_grid->Mx() != icebin_specI->nx()) || (pism_grid->My() != icebin_specI->ny())) {
            (*icebin_error)(-1,
                "Grid mismatch: pism=(%d, %d) icebin=(%d, %d)", pism_grid->Mx(), pism_grid->My(), icebin_specI->nx(), icebin_specI->ny());
        }
    }

    printf("END IceCoupler_PISM::_cold_start()\n");
}

void IceCoupler_PISM::run_timestep(double time_s,
    blitz::Array<double,2> const &ice_ivalsI,    // ice_ivalsI(nvar, nI)
    blitz::Array<double,2> &ice_ovalsI,    // ice_ovalsI(nvar, nI)
    bool run_ice)    // Should we run the ice model?
{
    PetscErrorCode ierr;

    printf("BEGIN IceCoupler_PISM::run_timestep()\n");

    // ----------- Bounds Checking
    // Check dimensions
    std::array<long,5> extents0{
        ice_ivalsI.extent(0),
        pism_ivars.size(),
        ice_ovalsI.extent(0),
        pism_ovars.size(),
        ice_ivalsI.extent(1)};
    std::array<long,5> extents1{
        contract[IceCoupler::INPUT].size(),
        contract[IceCoupler::INPUT].size(),
        contract[IceCoupler::OUTPUT].size(),
        contract[IceCoupler::OUTPUT].size(),
        ice_ovalsI.extent(1)};
    if (extents0 != extents1) (*icebin_error)(-1,
        "Extents mismatch (%d=%d), (%d=%d), (%d=%d), (%d=%d), (%d=%d)",
        extents0[0], extents1[0],
        extents0[1], extents1[1],
        extents0[2], extents1[2],
        extents0[3], extents1[3],
        extents0[4], extents1[4]);
    size_t const nI = ice_ivalsI.extent(1);

    // Check Petsc types
    if (sizeof(double) != sizeof(PetscScalar)) {
        (*icebin_error)(-1, "PetscScalar must be same as double\n");
    }

    if (run_ice) {
        // ---------- Load input into PISM's PETSc arrays
        // Fill pism_ivars[i] <-- iceIvals[:,i]
        // pism_ivars are distributed (global) vectors.
        for (int ivar=0; ivar<contract[IceCoupler::INPUT].size(); ++ivar) {
            VarMeta const &cf(contract[IceCoupler::INPUT][ivar]);

            // -------- Use
            // NOTE: vtmp_p0_va and *vtmp_p0 have different types
            pism::petsc::VecArray vtmp_p0_va(*vtmp_p0);
            if (am_i_root()) {
                blitz::Array<double,1> vtmp_b(vtmp_p0_va.get(), blitz::shape(nI), blitz::neverDeleteData);
                for (int iI=0; iI<nI; ++iI) {
                    vtmp_b(iI) = ice_ivalsI(ivar, iI);
                }
            }
            if (!(cf.flags & contracts::PRIVATE)) {
                IceModelVec2S *pism_var = pism_ivars[ivar];
                pism_var->get_from_proc0(*vtmp_p0);
            }
        }

        // -------- Figure out the timestep
        pism_in_nc->write(time_s);


        // =========== Run PISM for one coupling timestep
        // Time of last time we coupled
        auto time0(pism_grid->ctx()->time()->current());
        printf("BEGIN pism_ice_model->run_to(%f -> %f) %p\n",
            time0, time_s, pism_ice_model.get());
        // See pism::icebin::IBIceModel::run_to()
        pism_ice_model->run_to(time_s);
        printf("END pism_ice_model->run_to()\n");
        auto time1(pism_grid->ctx()->time()->current());


        if ((pism_ice_model->mass_t() != time_s) || (pism_ice_model->enthalpy_t() != time_s) || (time1 != time_s)) {
            (*icebin_error)(-1,
                "ERROR: PISM time (mass=%f, enthalpy=%f) doesn't match ICEBIN time %f",
                pism_ice_model->mass_t(), pism_ice_model->enthalpy_t(), time_s);
        }

        pism_ice_model->set_rate(time1 - time0);
    }    // if run_ice

    pism_ice_model->prepare_outputs(time_s);

    pism_out_nc->write(time_s);    // Writes from PISM-format variables on I grid

    get_state(ice_ovalsI, run_ice ? 0 : contracts::INITIAL);

    pism_ice_model->reset_rate();
    printf("END IceCoupler_PISM::run_timestep()\n");
}

/** Copies PISM->Icebin output variables from PISM variables to
the Icebin-supplied variables (on the root node).
@param mask Only do it for variables where (flags & mask) == mask.  Set to 0 for "all." */
void IceCoupler_PISM::get_state(
    blitz::Array<double,2> &ice_ovalsI,    // ice_ovalsI(nvar, nI)
    unsigned int mask)
{
    printf("BEGIN IceCoupler_PISM::get_state: %ld (mask = %d)\n", pism_ovars.size(), mask);
    VarSet const &ocontract(contract[IceCoupler::OUTPUT]);

    // Copy the outputs to the blitz arrays
    int nI = ice_ovalsI.extent(1);
    for (unsigned int ivar=0; ivar<pism_ovars.size(); ++ivar) {
        blitz::Array<double,1> ice_ovalsI_ivar(nI);
        VarMeta const &cf(ocontract.data[ivar]);
        bool const priv = (cf.flags & contracts::PRIVATE);

        if (!(pism_ovars[ivar] || priv)) (*icebin_error)(-1,
            "IceCoupler_PISM: Contract output %s (modele_pism.cpp) is not linked up to a pism_ovar (MassEnergyBudget.cpp)", ocontract.index[ivar].c_str());

        if ((cf.flags & mask) != mask) continue;

        printf("IceCoupler_PISM::get_state(mask=%d) copying field %s\n", mask, cf.name.c_str());

        if (am_i_root()) {      // ROOT in PISM communicator
            // Get matching input (val) and output (pism_var) variables
            iceModelVec2S_to_blitz_xy(*pism_ovars[ivar], ice_ovalsI_ivar);    // Allocates oval2_xy if needed

            // Copy to the output array
            ice_ovalsI(ivar, blitz::Range::all()) = ice_ovalsI_ivar;
        } else {
            blitz::Array<double,2> oval2_xy;    // dummy
            iceModelVec2S_to_blitz_xy(*pism_ovars[ivar], ice_ovalsI_ivar);
        }

        // Now send those data from the PISM root to the GCM root (MPI nodes)
        // (DUMMY for now, just make sure PISM and GCM have the same root)
        if (pism_root != gcm_coupler->gcm_params.gcm_root) (*icebin_error)(-1,
            "PISM and the GCM must share the same root!");
    }
    printf("END IceCoupler_PISM::get_state\n");
}


// ============================================================================
// Utility Functions...

void IceCoupler_PISM::deallocate()
{
    PetscErrorCode ierr;

//    ierr = VecDestroy(&g2); PISM_CHK(ierr, "VecDestroy");
//    ierr = VecDestroy(&g2natural); PISM_CHK(ierr, "VecDestroy");
    // ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);
//  ierr = VecDestroy(&Hp0); CHKERRQ(ierr);
}


void IceCoupler_PISM::iceModelVec2S_to_blitz_xy(
    pism::IceModelVec2S const &pism_var,
    blitz::Array<double,1> &ret1)    // 1D version of 2D array
{
    PetscErrorCode ierr;

    if (am_i_root()) {
        auto xy_shape(blitz::shape(ny(), nx()));
        if (ret1.size() == 0) {
            ret1.reference(blitz::Array<double,1>(ny() * nx()));
        } else {
            if (ret1.extent(0) != xy_shape[0] * xy_shape[1])
                (*icebin_error)(-1,
                "IceCoupler_PISM::iceModelVec2S_to_blitz_xy(): "
                "ret1(%d) should be (%d * %d)",
                ret1.extent(0), xy_shape[0], xy_shape[1]);
        }
    }

#if 0
// Don't know how to translate get_dof() to current PISM API
    if (pism_var.get_dof() != 1)
        SETERRQ(pism_grid->com, 1, "This method only supports IceModelVecs with dof == 1");
#endif

    // Gather data to one processor
    // OLD CODE: PetscScalar **bHp0;

    // OLD CODE: pism_var.put_on_proc0(Hp0, scatter, g2, g2natural);
    // 
    // See ec80c25e (2014-10-08)
    // Rewrite put_on_proc0/get_from_proc0 IceModelVec2S methods.
    // 
    // Also add IceModelVec2S::allocate_proc0_copy(), which allocates scatters
    // and work space the first time it is called for a given DM and re-uses
    // them for all IceModelVec2Ss using this DM.
    // 
    // To transfer a 2D field to processor 0 and back, do this:
    // 
    // IceModelVec2S field;
    // Vec field_p0;
    // 
    // field.create(...);
    // field.allocate_proc0_copy(field_p0);
    // field.put_on_proc0(field_p0);
    // // do work on processor 0
    // field.get_from_proc0(field_p0);
    // VecDestroy(&field_p0);

    Hp0 = pism_var.allocate_proc0_copy();
    pism_var.put_on_proc0(*Hp0);


    // Copy it to blitz array (on the root node only)
    // OLD CODE: ierr = VecGetArray2d(Hp0, pism_grid->Mx, pism_grid->My, 0, 0, &bHp0);

//    if (am_i_root()) ret1 = nan;  // Vector operation, initializes ice_ovalsI

//    PetscInt Hp0_size;
//    VecGetLocalSize(*Hp0, &Hp0_size);



//    petsc::VecArray va(*Hp0);
//    double const *Hp0_data = va.get();
    if (am_i_root()) {
        double const *Hp0_data = petsc::VecArray(*Hp0).get();

        long nI = pism_grid->Mx() * pism_grid->My();
        ret1 = nan;
        for (int i=0; i<nI; ++i) ret1(i) = Hp0_data[i];
    }
}


void IceCoupler_PISM::transfer_constant(std::string const &dest, std::string const &src, double multiply_by, bool set_new)
{
    // Make sure the PISM constant already exists
    if (!set_new && !pism_config()->is_set(dest)) (*icebin_error)(-1,
        "IceCoupler_PISM::transfer_constant: Trying to set '%s', which is not a PISM configuration parameter.  Is it misspelled?", dest.c_str());

    // Discover the units PISM requires.
    std::string units = pism_config()->get_string(dest + "_units");
    double val = gcm_coupler->gcm_constants.get_as(src, units) * multiply_by;
    pism_config()->set_double(dest, val);
printf("IceCoupler_PISM::transfer_constant: %s = %g [%s] (from %s in GCM)\n", dest.c_str(), val, units.c_str(), src.c_str());




}

void IceCoupler_PISM::set_constant(std::string const &dest, double src_val, std::string const &src_units, bool set_new)
{
    // Make sure the PISM constant already exists
    if (!set_new && !pism_config()->is_set(dest)) (*icebin_error)(-1,
        "IceCoupler_PISM::set_constant: Trying to set '%s', which is not a PISM configuration parameter.  Is it misspelled?", dest.c_str());

    ibmisc::ConstantSet const &gcm_constants(gcm_coupler->gcm_constants);

    // Discover the units PISM requires.
    std::string dest_units = pism_config()->get_string(dest + "_units");

    UTUnit usrc(gcm_constants.ut_system->parse(src_units));
    UTUnit udest(gcm_constants.ut_system->parse(dest_units));
    CVConverter cv(usrc, udest);
    double dest_val = cv.convert(src_val);

    pism_config()->set_double(dest, dest_val);
printf("IceCoupler_PISM::transfer_constant: %s = %g %s (from %s in GCM)\n", dest.c_str(), dest_val, dest_units.c_str(), usrc.c_str());
}


}}
