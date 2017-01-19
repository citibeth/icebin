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

void IceCoupler_PISM::ncread(ibmisc::NcIO &ncio, std::string const &vname_sheet)
{
    IceCoupler::ncread(ncio, vname_sheet);

    printf("BEGIN IceCoupler_PISM::ncread()\n");
    GCMParams const &_gcm_params(gcm_coupler->gcm_params);

    this->icebin_gridI = dynamic_cast<Grid_XY const *>(gridI());

    // General args passed to the ice sheet, regardless of which ice model is being used
    NcVar info_var(ncio.nc->getVar(vname_sheet + ".info"));
    // PISM parameters, passed to PISM via argv
    NcVar pism_var(ncio.nc->getVar(vname_sheet + ".pism"));

    // Get simple arguments
    ibmisc::get_or_put_att(info_var, 'r', "update_elevation", &update_elevation, 1);

    // Create arguments from PISM configuration
    pism_args.push_back("icebin_pism");

    // Get arguments from IceBin configuration
    std::map<std::string, NcVarAtt> pism_atts(pism_var.getAtts());
    for (auto jj=pism_atts.begin(); jj != pism_atts.end(); ++jj) {
        std::string const &name(jj->first);
        NcAtt const &att(jj->second);
        std::string val;
        att.getValues(val);

        auto ii(path_args.find(name));
        if (ii == path_args.end()) {
            // Regular case, just use the value there
        } else {
            // Resolve path names according to the configuration directory
            auto &resolve(ii->second == "i" ?
                gcm_coupler->gcm_params.config_dir : gcm_coupler->gcm_params.run_dir);

            val = boost::filesystem::absolute(
                boost::filesystem::path(val),
                resolve).string();
            printf("IceCoupler_PISM resolving %s: --> %s\n", name.c_str(), val.c_str());
        }

        pism_args.push_back("-" + name);
        pism_args.push_back(val);
    }

#if 0
    // Hard-code these variables because we will always need them
    pism_args.push_back("-extra_vars");
    pism_args.push_back("climatic_mass_balance_cumulative,nonneg_flux_cumulative,grounded_basal_flux_cumulative,floating_basal_flux_cumulative,flux_divergence");
#endif
}
// ======================================================================
void IceCoupler_PISM::cold_start(
    ibmisc::Datetime const &time_base,
    double time_start_s)
{
    // Call overridden method
    IceCoupler::cold_start(time_base, time_start_s);

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
for (int i=0; i<argc; ++i) printf(" %s", argv[i]);
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

    // Initialize Petsc
    // petsc_initializer.reset(new pism::petsc::Initializer(pism_comm, argc, argv));
    petsc_initializer.reset(new pism::petsc::Initializer(argc, argv, "IceBin GCM Coupler"));

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
    config->set_string("reference_date", reference_date);
    // ------------------------------ //

#if 0    // Don't bother with profiling inside of GCM
    if (profiling_log.is_set()) {
      ctx->profiling().start();
    }
#endif

    log->message(3, "* Setting the computational grid...\n");
    IceGrid::Ptr pism_grid = IceGrid::FromOptions(ctx);

    pism::icebin::IBIceModel::Params params;
        params.time_start_s = gcm_coupler->time_start_s;
        params.output_dir = boost::filesystem::absolute(gcm_coupler->gcm_params.run_dir).string();
    pism_ice_model.reset(new pism::icebin::IBIceModel(pism_grid, ctx, params));

    // ------------------------------------------- \\

    // Transfer constants from GCM to PISM, and also set up coupling contracts.
    // This is the right place to do it, since the PISM systme is fully up and functional,
    // and all PISM config files have been read.
    // This call through the GCMCoupler will call back to setup_contracts_xxx().
    contracts::setup(*gcm_coupler, *this);

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
        pism_ivars[ix] = &pism_surface_model->icebin_massxfer;

    // Ignore surface_temp, it is not useful...
    ix = contract[INPUT].index.at("enthxfer");
        pism_ivars[ix] = &pism_surface_model->icebin_enthxfer;

    ix = contract[INPUT].index.at("deltah");
        pism_ivars[ix] = &pism_surface_model->icebin_deltah;

    // Check that all PISM inputs are bound to a variable
    bool err = false;
    for (unsigned int i=0; i<pism_ivars.size(); ++i) {
        IceModelVec2S *pism_var = pism_ivars[i];
        if (!pism_var) fprintf(stderr,
            "PISM input %s is not bound to a variable\n",
            contract[INPUT].data[ix].name.c_str());
    }
    if (err) (*icebin_error)(-1, "Exiting due to errors");


    // Initialize scatter/gather stuff
    da2 = pism_grid->get_dm(1, // dm_dof
        pism_grid->ctx()->config()->get_double("grid_max_stencil_width"));

//    ierr = DMCreateGlobalVector(*da2, &g2); PISM_CHK(ierr, "DMCreateGlobalVector");

    // note we want a global Vec but reordered in the natural ordering
    // so when it is scattered to proc zero it is not all messed up;
    // see above
//    ierr = DMDACreateNaturalVector(*da2, &g2natural);
//        PISM_CHK(ierr, "DMDACreateNaturalVector");

    // next get context *and* allocate samplep0 (on proc zero only, naturally)
//    ierr = VecScatterCreateToZero(g2natural, &scatter, &Hp0);
//        PISM_CHK(ierr, "VecScatterCreateToZero");


    // ============== Set up variables for OUTPUT contract
    // -------------- Allocate blitz:Array<double,1> output variables,
    // which are passed back to Icebin.

    // -------------- Link to PISM-format output variables, used to fill ovars
    pism_ovars.resize(contract[OUTPUT].size(), NULL);
    ix = contract[OUTPUT].index.at("ice_surface_elevation");       // Elevation of top surface of ice sheet
        pism_ovars[ix] = &pism_ice_model->ice_surface_elevation(); // see PISM's iceModel.hh

    ix = contract[OUTPUT].index.at("ice_surface_elevation");
        pism_ovars[ix] = &pism_ice_model->ice_surface_elevation();
    ix = contract[OUTPUT].index.at("ice_thickness");
        pism_ovars[ix] = &pism_ice_model->ice_thickness();
    ix = contract[OUTPUT].index.at("bed_topography");
        pism_ovars[ix] = &pism_ice_model->bed_model()->bed_elevation();

    ix = contract[OUTPUT].index.at("mask");
        pism_ovars[ix] = &pism_ice_model->cell_type();

    // Mass of top two layers
    ix = contract[OUTPUT].index.at("M1");
        pism_ovars[ix] = &pism_ice_model->M1;
    ix = contract[OUTPUT].index.at("M2");
        pism_ovars[ix] = &pism_ice_model->M2;

    // Enthalpy of top two layers
    ix = contract[OUTPUT].index.at("H1");
        pism_ovars[ix] = &pism_ice_model->H1;
    ix = contract[OUTPUT].index.at("H2");
        pism_ovars[ix] = &pism_ice_model->H2;

    // Volume of top two layers
    ix = contract[OUTPUT].index.at("V1");
        pism_ovars[ix] = &pism_ice_model->V1;
    ix = contract[OUTPUT].index.at("V2");
        pism_ovars[ix] = &pism_ice_model->V2;

    // For MassEnergyBudget variables that have a contract name specified,
    // link them up into pism_ovars now.
    for (auto ii = pism_ice_model->rate.all_vecs.begin(); ii != pism_ice_model->rate.all_vecs.end(); ++ii) {
        if (ii->contract_name == "") continue;

        int ix = contract[OUTPUT].index.at(ii->contract_name);
        pism_ovars[ix] = &ii->vec;
    }

    // -------------- Initialize pism_out.nc
    {
        boost::filesystem::path output_dir(params.output_dir);
        std::string ofname = (output_dir / "pism_out.nc").string();
        std::vector<pism::IceModelVec const *> vecs;
        for (auto &vec : pism_ovars) vecs.push_back(vec);
        pism_out_nc.reset(new pism::icebin::VecBundleWriter(
            pism_ice_model->grid(), ofname, vecs));
        pism_out_nc->init();
    }

    // ------------- Initialize pism_in.nc
    {
        boost::filesystem::path output_dir(params.output_dir);
        std::string ofname = (output_dir / "pism_in.nc").string();
        std::vector<pism::IceModelVec const *> vecs;
        for (auto &vec : pism_ivars) vecs.push_back(vec);
        pism_in_nc.reset(new pism::icebin::VecBundleWriter(
            pism_ice_model->grid(), ofname, vecs));
        pism_in_nc->init();
    }


    // ============== Miscellaneous
    // Check that grid dimensions match
    if ((pism_grid->Mx() != icebin_gridI->nx()) || (pism_grid->My() != icebin_gridI->ny())) {
        (*icebin_error)(-1,
            "Grid mismatch: pism=(%d, %d) icebin=(%d, %d)", pism_grid->Mx(), pism_grid->My(), icebin_gridI->nx(), icebin_gridI->ny());
    }

    printf("END IceCoupler_PISM::allocate()\n");
}

blitz::Array<double,1> IceCoupler_PISM::get_elevI()
{
    // Reshape 1D blitz variable to 2D for use with PISM
    blitz::Array<double,1> retI(ny()*nx());
    auto ret_xy(reshape<double,1,2>(retI, blitz::shape(ny(), nx())));

    iceModelVec2S_to_blitz_xy(
        pism_ice_model->ice_surface_elevation(),
        ret_xy);
    return retI;
}

void IceCoupler_PISM::run_timestep(double time_s,
    blitz::Array<double,2> const &ice_ivalsI,    // ice_ivalsI(nI, nvar)
    blitz::Array<double,2> const &ice_ovalsI,    // ice_ovalsI(nI, nvar)
    bool run_ice,    // Should we run the ice model?
    bool am_i_root)
{
    PetscErrorCode ierr;

    // ----------- Type Checking
    // Check dimensions
    std::array<long,5> extents0{
        ice_ivalsI.extent(1),
        pism_ivars.size(),
        ice_ovalsI.extent(1),
        pism_ovars.size(),
        ice_ivalsI.extent(0)};
    std::array<long,5> extents1{
        contract[IceCoupler::INPUT].size(),
        contract[IceCoupler::INPUT].size(),
        contract[IceCoupler::OUTPUT].size(),
        contract[IceCoupler::OUTPUT].size(),
        ice_ovalsI.extent(0)};
    if (extents0 != extents1) (*icebin_error)(-1,
        "Extents mismatch (%d=%d), (%d=%d), (%d=%d), (%d=%d), (%d=%d)",
        extents0[0], extents1[0],
        extents0[1], extents1[1],
        extents0[2], extents1[2],
        extents0[3], extents1[3],
        extents0[4], extents1[4]);

    // Check Petsc types
    if (sizeof(double) != sizeof(PetscScalar)) {
        (*icebin_error)(-1, "PetscScalar must be same as double\n");
    }

    size_t nI = ice_ivalsI.extent(0);

    if (run_ice) {
        // ---------- Load input into PISM's PETSc arrays
        // Fill pism_ivars[i] <-- iceIvals[:,i]
        // pism_ivars are distributed (global) vectors.
        blitz::Array<PetscScalar,1> g2_y(nI);
        blitz::Array<int,1> g2_ix(nI);
        for (int ivar=0; ivar<contract[IceCoupler::INPUT].size(); ++ivar) {
            VarMeta const &cf(contract[IceCoupler::INPUT][ivar]);

            // Get matching input (val) and output (pism_var) variables
            IceModelVec2S *pism_var = pism_ivars[ivar];

            // Inputs specified in the contract are not (necessarily) attached
            // to any PISM var.  If they are not, just drop them on the ground.
            if (!pism_var) continue;

            // Copy value to a stride=1 array
            for (int iI=0; iI<nI; ++iI) {
                g2_ix(iI) = iI;
                g2_y(iI) = ice_ivalsI(iI, ivar);
            }

            // Put into a natural-ordering global distributed Petsc Vec
            ierr = VecSet(g2natural, cf.default_value); PISM_CHK(ierr, "run_timestep");
            ierr = VecSetValues(g2natural, nI, &g2_ix(0), &g2_y(0), INSERT_VALUES); PISM_CHK(ierr, "run_timestep");

            ierr = VecAssemblyBegin(g2natural); PISM_CHK(ierr, "run_timestep");
            ierr = VecAssemblyEnd(g2natural); PISM_CHK(ierr, "run_timestep");

            // Copy to Petsc-ordered global vec
            ierr = DMDANaturalToGlobalBegin(*da2, g2natural, INSERT_VALUES, g2); PISM_CHK(ierr, "run_timestep");
            ierr = DMDANaturalToGlobalEnd(*da2, g2natural, INSERT_VALUES, g2); PISM_CHK(ierr, "run_timestep");

            // Copy to the output variable
            // (Could we just do DMDANaturalToGlobal() directly to this?)
            pism_var->copy_from_vec(g2);
        }

        // -------- Figure out the timestep
        auto old_pism_time(pism_grid->ctx()->time()->current()); // beginning of this PISM timestep [s]
        auto timestep_s = time_s - old_pism_time;       // [s]

        // -------- Determine Dirichlet B.C. for ice sheet
        // This is done by taking the changes in the "borrowed" enthalpies
        // from the GCM, and applying them to the corresponding top layer in
        // the ice model.  The result is placed in surface->surface_temp.
        pism::icebin::IBSurfaceModel * const surface = pism_ice_model->ib_surface_model();
        pism_ice_model->construct_surface_temp(
            surface->icebin_deltah,
            contract[INPUT].at("deltah").default_value,
            timestep_s,
            surface->surface_temp);

        pism_in_nc->write(time_s);


        // =========== Run PISM for one coupling timestep
        // Time of last time we coupled
        printf("BEGIN pism_ice_model->run_to(%f -> %f) %p\n",
            pism_grid->ctx()->time()->current(), time_s, pism_ice_model.get());
        // See pism::icebin::IBIceModel::run_to()
        pism_ice_model->run_to(time_s);
        printf("END pism_ice_model->run_to()\n");


        if ((pism_ice_model->mass_t() != time_s) || (pism_ice_model->enthalpy_t() != time_s)) {
            (*icebin_error)(-1,
                "ERROR: PISM time (mass=%f, enthalpy=%f) doesn't match ICEBIN time %f", pism_ice_model->mass_t(), pism_ice_model->enthalpy_t(), time_s);
        }

        pism_ice_model->set_rate(pism_ice_model->enthalpy_t() - old_pism_time);
    }    // if run_ice

    pism_ice_model->prepare_outputs(time_s);

    pism_out_nc->write(time_s);

    get_state(ice_ovalsI, run_ice ? contracts::INITIAL : 0);
    pism_ice_model->reset_rate();
}

/** Copies PISM->Icebin output variables from PISM variables to
the Icebin-supplied variables (on the root node).
@param mask Only do it for variables where (flags & mask) == mask.  Set to 0 for "all." */
void IceCoupler_PISM::get_state(
    blitz::Array<double,2> const &ice_ovalsI,    // ice_ovalsI(nI, nvar)
    unsigned int mask)
{
    printf("BEGIN IceCoupler_PISM::get_state: %ld\n", pism_ovars.size());
    VarSet const &ocontract(contract[IceCoupler::OUTPUT]);

    // Copy the outputs to the blitz arrays
    for (unsigned int ivar=0; ivar<pism_ovars.size(); ++ivar) {
        blitz::Array<double,1> const &ice_ovalsI_ivar(
            ice_ovalsI(blitz::Range::all(),ivar));

        if (!pism_ovars[ivar]) (*icebin_error)(-1,
            "IceCoupler_PISM: Contract output %s (modele_pism.cpp) is not linked up to a pism_ovar (MassEnergyBudget.cpp)", ocontract.index[ivar].c_str());

        VarMeta const &cf(ocontract.data[ivar]);
        if ((cf.flags & mask) != mask) continue;

        printf("IceCoupler_PISM::get_state(mask=%d) copying field %s\n", mask, cf.name.c_str());

        if (am_i_root()) {      // ROOT in PISM communicator

            // Reshape 1D blitz variable to 2D for use with PISM
            blitz::Array<double,2> oval2_xy(
                reshape<double,1,2>(ice_ovalsI_ivar, blitz::shape(ny(), nx())));
        
    
            // Get matching input (val) and output (pism_var) variables
            iceModelVec2S_to_blitz_xy(*pism_ovars[ivar], oval2_xy);    // Allocates oval2_xy if needed

        } else {
            blitz::Array<double,2> oval2_xy;    // dummy
            iceModelVec2S_to_blitz_xy(*pism_ovars[ivar], oval2_xy);
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
    blitz::Array<double,2> &ret)
{
    PetscErrorCode ierr;

    if (am_i_root()) {
        auto xy_shape(blitz::shape(nx(), ny()));
        if (ret.size() == 0) {
            ret.reference(blitz::Array<double,2>(ny(), nx()));
        } else {
            if (ret.extent(0) != xy_shape[0] || ret.extent(1) != xy_shape[1])
                (*icebin_error)(-1,
                "IceCoupler_PISM::iceModelVec2S_to_blitz_xy(): "
                "ret(%d, %d) should be (%d, %d)",
                ret.extent(0), ret.extent(1), xy_shape[0], xy_shape[1]);
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

    if (am_i_root()) {
        // Copy it to blitz array (on the root node only)
        // OLD CODE: ierr = VecGetArray2d(Hp0, pism_grid->Mx, pism_grid->My, 0, 0, &bHp0);

        ret = nan;  // Vector operation, initializes ice_ovalsI
        for (PetscInt i=0; i < pism_grid->Mx(); i++) {
            for (PetscInt j=0; j < pism_grid->My(); j++) {
// TODO: ???? What should this be?
//                ret(i, j) = bHp0[i][j];
            }
        }
        // OLD CODE: ierr = VecRestoreArray2d(Hp0, pism_grid->Mx, pism_grid->My, 0, 0, &bHp0); CHKERRQ(ierr);
    }
}


void IceCoupler_PISM::transfer_constant(std::string const &dest, std::string const &src, double multiply_by, bool set_new)
{
    // Make sure the PISM constant already exists
    if (!set_new && !pism_config()->is_set(dest)) (*icebin_error)(-1,
        "IceCoupler_PISM::transfer_constant: Trying to set '%s', which is not a PISM configuration parameter.  Is it misspelled?", dest.c_str());

    // Discover the units PISM requires.
    std::string doc = pism_config()->get_string(dest + "_doc");
    std::string units = doc.substr(0, doc.find(';'));
    double val = gcm_coupler->gcm_constants.get_as(src, units) * multiply_by;
    pism_config()->set_double(dest, val);
printf("IceCoupler_PISM::transfer_constant: %s = %g %s (from %s in GCM)\n", dest.c_str(), val, units.c_str(), src.c_str());
}

void IceCoupler_PISM::set_constant(std::string const &dest, double src_val, std::string const &src_units, bool set_new)
{
    // Make sure the PISM constant already exists
    if (!set_new && !pism_config()->is_set(dest)) (*icebin_error)(-1,
        "IceCoupler_PISM::set_constant: Trying to set '%s', which is not a PISM configuration parameter.  Is it misspelled?", dest.c_str());

    ibmisc::ConstantSet const &gcm_constants(gcm_coupler->gcm_constants);

    // Discover the units PISM requires.
    std::string doc = pism_config()->get_string(dest + "_doc");
    std::string dest_units = doc.substr(0, doc.find(';'));

    UTUnit usrc(gcm_constants.ut_system->parse(src_units));
    UTUnit udest(gcm_constants.ut_system->parse(dest_units));
    CVConverter cv(usrc, udest);
    double dest_val = cv.convert(src_val);

    pism_config()->set_double(dest, dest_val);
printf("IceCoupler_PISM::transfer_constant: %s = %g %s (from %s in GCM)\n", dest.c_str(), dest_val, dest_units.c_str(), usrc.c_str());
}


}}
