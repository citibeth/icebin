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

#include <base/stressbalance/PISMStressBalance.hh>
#include <earth/PISMBedDef.hh>

#include <spsparse/sort.hpp>

#include <ibmisc/netcdf.hpp>
#include <ibmisc/ibmisc.hpp>

#include <icebin/GCMCoupler.hpp>
#include <icebin/pism/IceModel_PISM.hpp>
#include <icebin/contracts/contracts.hpp>

extern "C" void libpismutil_refaddr();

using namespace ibmisc;
using namespace pism;
using namespace netCDF;

namespace icebin {
namespace gpism {

static double const nan = std::numeric_limits<double>::quiet_NaN();

IceModel_PISM::IceModel_PISM()
    : IceModel(IceModel::Type::PISM),
    write_pism_inputs(true)
{
    printf("BEGIN/END IceModel_PISM::IceModel_PISM\n");
}


void IceModel_PISM::update_ice_sheet(ibmisc::NcIO &ncio, std::string const &vname)
{
    printf("BEGIN IceModel_PISM::update_ice_sheet(%s)\n", vname.c_str());

    auto pism_var(ncio.nc->getVar(vname + ".pism"));    // PISM Parameters
    std::string pism_i_str;
    get_or_put_att(pism_var, 'r', "i", pism_i_str, true);  // PISM -i argument (input file)

    std::string pism_i = boost::filesystem::absolute(
        boost::filesystem::path(pism_i_str),
        coupler->gcm_params.config_dir).string();

    // Read variables from PISM input file
    // byte mask(time, x, y) ;
    //      mask:units = "" ;
    //      mask:coordinates = "lat lon" ;
    //      mask:flag_meanings = "ice_free_bedrock grounded_ice floating_ice ice_free_ocean" ;
    //      mask:grid_mapping = "mapping" ;
    //      mask:long_name = "ice-type (ice-free/grounded/floating/ocean) integer mask" ;
    //      mask:pism_intent = "diagnostic" ;
    //      mask:flag_values = 0b, 2b, 3b, 4b ;
    // double thk(time, x, y) ;
    //      thk:units = "m" ;
    //      thk:valid_min = 0. ;
    //      thk:coordinates = "lat lon" ;
    //      thk:grid_mapping = "mapping" ;
    //      thk:long_name = "land ice thickness" ;
    //      thk:pism_intent = "model_state" ;
    //      thk:standard_name = "land_ice_thickness" ;
    // double topg(time, x, y) ;
    //      topg:units = "m" ;
    //      topg:coordinates = "lat lon" ;
    //      topg:grid_mapping = "mapping" ;
    //      topg:long_name = "bedrock surface elevation" ;
    //      topg:pism_intent = "model_state" ;
    //      topg:standard_name = "bedrock_altitude" ;
printf("Opening PISM file for elev2 and mask2: %s\n", pism_i.c_str());
    NcIO ncin(pism_i, netCDF::NcFile::read);
    size_t ntime = ncin.nc->getDim("time").getSize();
    size_t nx = ncin.nc->getDim("x").getSize();
    size_t ny = ncin.nc->getDim("y").getSize();

    blitz::Array<unsigned char,2> mask(nx, ny);
    NcVar mask_var = ncin.nc->getVar("mask");
    mask_var.getVar({ntime-1, 0, 0}, {1, nx, ny}, mask.data());

    blitz::Array<double,2> thk(nx, ny);
    NcVar thk_var = ncin.nc->getVar("thk");
    thk_var.getVar({ntime-1, 0, 0}, {1, nx, ny}, thk.data());

    blitz::Array<double,2> topg(nx, ny);
    NcVar topg_var = ncin.nc->getVar("topg");
    topg_var.getVar({ntime-1, 0, 0}, {1, nx, ny}, topg.data());

    ncin.close();

    // Set update_elevation=false to temporarily "repair" fields that are
    // broken due to a difference in the elevation classes used to generate
    // the SMB in a one-way coupled run, and the elevation classes
    // as defined by the ice sheet.
    if (update_elevation) {
        sheet->elevI.clear();
        for (int i=0; i<nx; ++i) {
        for (int j=0; j<ny; ++j) {
            if (mask(i,j) == 2) {
                long ix2 = gridI()->indexing.tuple_to_index<3>({i,j});
                double elev = topg(i,j) + thk(i,j);
                sheet->elevI.add({ix2}, elev);
            }
        }}
    }
    printf("END IceModel_PISM::update_ice_sheet()\n");
}



// See: http://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

void IceModel_PISM::transfer_constant(std::string const &dest, std::string const &src, double multiply_by, bool set_new)
{
    auto &config(*ice_model->ctx()->config());   // PISM configuration

    // Make sure the PISM constant already exists
    if (!set_new && !config.is_set(dest)) (*icebin_error)(-1,
        "IceModel_PISM::transfer_constant: Trying to set '%s', which is not a PISM configuration parameter.  Is it misspelled?", dest.c_str());

    // Discover the units PISM requires.
    std::string doc = config.get_string(dest + "_doc");
    std::string units = doc.substr(0, doc.find(';'));
    double val = coupler->gcm_constants.get_as(src, units) * multiply_by;
    config.set_double(dest, val);
printf("IceModel_PISM::transfer_constant: %s = %g %s (from %s in GCM)\n", dest.c_str(), val, units.c_str(), src.c_str());
}

void IceModel_PISM::set_constant(std::string const &dest, double src_val, std::string const &src_units, bool set_new)
{
    auto &config(*ice_model->ctx()->config());   // PISM configuration

    // Make sure the PISM constant already exists
    if (!set_new && !config.is_set(dest)) (*icebin_error)(-1,
        "IceModel_PISM::set_constant: Trying to set '%s', which is not a PISM configuration parameter.  Is it misspelled?", dest.c_str());

    ibmisc::ConstantSet const &gcm_constants(coupler->gcm_constants);

    // Discover the units PISM requires.
    std::string doc = config.get_string(dest + "_doc");
    std::string dest_units = doc.substr(0, doc.find(';'));

    UTUnit usrc(gcm_constants.ut_system->parse(src_units));
    UTUnit udest(gcm_constants.ut_system->parse(dest_units));
    CVConverter cv(usrc, udest);
    double dest_val = cv.convert(src_val);

    config.set_double(dest, dest_val);
printf("IceModel_PISM::transfer_constant: %s = %g %s (from %s in GCM)\n", dest.c_str(), dest_val, dest_units.c_str(), usrc.c_str());
}



// Arguments that are paths, and thus need pathname resolution
// For stable0.5 branch
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

void IceModel_PISM::ncread(ibmisc::NcIO &ncio, std::string const &vname_sheet)
{
    IceModel::ncread(ncio, vname_sheet);

    printf("BEGIN IceModel_PISM::ncread()\n");
    GCMParams const &_gcm_params(coupler->gcm_params);

    this->icebin_gridI = dynamic_cast<Grid_XY const *>(gridI());

    // General args passed to the ice sheet, regardless of which ice model is being used
    NcVar info_var(ncio.nc->getVar(vname_sheet + ".info"));
    // PISM parameters, passed to PISM via argv
    NcVar pism_var(ncio.nc->getVar(vname_sheet + ".pism"));

    // Get simple arguments
    ibmisc::get_or_put_att(info_var, 'r', "update_elevtation", &update_elevation, 1);

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
                coupler->gcm_params.config_dir : coupler->gcm_params.run_dir);

            val = boost::filesystem::absolute(
                boost::filesystem::path(val),
                resolve).string();
            printf("IceModel_PISM resolving %s: --> %s\n", name.c_str(), val.c_str());
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

// -------------------------------------------

// Called from start_time_set().
void IceModel_PISM::allocate()
{
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
    pism_comm = coupler->gcm_params.gcm_comm;
    PetscErrorCode ierr;
    ierr = MPI_Comm_rank(pism_comm, &pism_rank); PISM_CHK(ierr, "MPI_Comm_rank");
    ierr = MPI_Comm_size(pism_comm, &pism_size); PISM_CHK(ierr, "MPI_Comm_size");

printf("[%d] pism_size = %d\n", pism_rank, pism_size);

    // Initialize Petsc
    // petsc_initializer.reset(new pism::petsc::Initializer(pism_comm, argc, argv));
    petsc_initializer.reset(new pism::petsc::Initializer(argc, argv, "IceBin GCM Coupler"));

    verbosityLevelFromOptions();
    Context::Ptr ctx = context_from_options(pism_comm, "IceModel_PISM");
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
    ibmisc::time::tm const &tb(coupler->gcm_params.time_base);
    std::string reference_date = (boost::format("%04d-%02d-%02d") % tb.year() % tb.month() % tb.mday()).str();
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
        params.time_start_s = coupler->gcm_params.time_start_s;
        params.output_dir = boost::filesystem::absolute(coupler->gcm_params.run_dir).string();
    ice_model.reset(new pism::icebin::IBIceModel(pism_grid, ctx, params));

    // ------------------------------------------- \\

    // Transfer constants from GCM to PISM, and also set up coupling contracts.
    // This is the right place to do it, since the PISM systme is fully up and functional,
    // and all PISM config files have been read.
    // This call through the GCMCoupler will call back to setup_contracts_xxx().
    contracts::setup(*coupler, *this);

//printf("start = %f\n", pism_grid->time->start());
//printf("end = %f\n", pism_grid->time->end());

    // This has the following stack trace:
    //  IceModel::init()                    [iceModel.cc]
    //  IceModel::model_state_setup()       [iMinit.cc]
    //  IceModel::init_couplers()           [iMinit.cc]
    //  surface->init()
    ice_model->init();

    // ============== Set up variables for INPUT contract

    // During the ice_model->init() call above the PISMIceModel
    // class (derived from PISM's IceModel) allocated an instance
    // of PSConstantICEBIN. This instance is owned and will be
    // de-allocated by PISMIceModel ice_model.
    pism_surface_model = ice_model->ib_surface_model();

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
        pism_ovars[ix] = &ice_model->ice_surface_elevation; // see PISM's iceModel.hh

    ix = contract[OUTPUT].index.at("ice_surface_elevation");
        pism_ovars[ix] = &ice_model->ice_surface_elevation;
    ix = contract[OUTPUT].index.at("ice_thickness");
        pism_ovars[ix] = &ice_model->ice_thickness;
    ix = contract[OUTPUT].index.at("bed_topography");
        pism_ovars[ix] = &ice_model->beddef->bed_elevation();

    ix = contract[OUTPUT].index.at("mask");
        pism_ovars[ix] = &ice_model->vMask;

    // Mass of top two layers
    ix = contract[OUTPUT].index.at("M1");
        pism_ovars[ix] = &ice_model->M1;
    ix = contract[OUTPUT].index.at("M2");
        pism_ovars[ix] = &ice_model->M2;

    // Enthalpy of top two layers
    ix = contract[OUTPUT].index.at("H1");
        pism_ovars[ix] = &ice_model->H1;
    ix = contract[OUTPUT].index.at("H2");
        pism_ovars[ix] = &ice_model->H2;

    // Volume of top two layers
    ix = contract[OUTPUT].index.at("V1");
        pism_ovars[ix] = &ice_model->V1;
    ix = contract[OUTPUT].index.at("V2");
        pism_ovars[ix] = &ice_model->V2;

    // For MassEnergyBudget variables that have a contract name specified,
    // link them up into pism_ovars now.
    for (auto ii = ice_model->rate.all_vecs.begin(); ii != ice_model->rate.all_vecs.end(); ++ii) {
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
            ice_model->grid(), ofname, vecs));
        pism_out_nc->init();
    }

    // ------------- Initialize pism_in.nc
    {
        boost::filesystem::path output_dir(params.output_dir);
        std::string ofname = (output_dir / "pism_in.nc").string();
        std::vector<pism::IceModelVec const *> vecs;
        for (auto &vec : pism_ivars) vecs.push_back(vec);
        pism_in_nc.reset(new pism::icebin::VecBundleWriter(
            ice_model->grid(), ofname, vecs));
        pism_in_nc->init();
    }


    // ============== Miscellaneous
    // Check that grid dimensions match
    if ((pism_grid->Mx() != icebin_gridI->nx()) || (pism_grid->My() != icebin_gridI->ny())) {
        (*icebin_error)(-1,
            "Grid mismatch: pism=(%d, %d) icebin=(%d, %d)", pism_grid->Mx(), pism_grid->My(), icebin_gridI->nx(), icebin_gridI->ny());
    }

    printf("END IceModel_PISM::allocate()\n");
}

void IceModel_PISM::deallocate()
{
    PetscErrorCode ierr;

//    ierr = VecDestroy(&g2); PISM_CHK(ierr, "VecDestroy");
//    ierr = VecDestroy(&g2natural); PISM_CHK(ierr, "VecDestroy");
    // ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);
//  ierr = VecDestroy(&Hp0); CHKERRQ(ierr);
}


// --------------------------------------------------------

// --------------------------------------------------
/** icebin_var Variable, already allocated, to receive data
@param icebin_var_xy The array to write into (on the root node).
If this array is not yet allocated (ROOT NODE ONLY), it will be allocated.*/
void IceModel_PISM::iceModelVec2S_to_blitz_xy(pism::IceModelVec2S const &pism_var,
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
                "IceModel_PISM::iceModelVec2S_to_blitz_xy(): "
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

        ret = nan;  // Vector operation, initializes ice_ovals_I
        for (PetscInt i=0; i < pism_grid->Mx(); i++) {
            for (PetscInt j=0; j < pism_grid->My(); j++) {
// TODO: ???? What should this be?
//                ret(i, j) = bHp0[i][j];
            }
        }
        // OLD CODE: ierr = VecRestoreArray2d(Hp0, pism_grid->Mx, pism_grid->My, 0, 0, &bHp0); CHKERRQ(ierr);
    }
}
// --------------------------------------------------
void IceModel_PISM::run_timestep(double time_s,
    blitz::Array<int,1> const &indices,
    std::vector<blitz::Array<double,1>> const &ivals2)
{
    PetscErrorCode ierr;

    printf("BEGIN IceModel_PISM::run_timestep(%f)\n", time_s);

    // This subroutine happens AFTER the gather to root.
    // THEREFORE... indices.size() will be 0, EXCEPT for root.

    int nvals = indices.extent(0);

    // Check number of variables matches for input
    if (ivals2.size() != pism_ivars.size()) (*icebin_error)(-1,
        "[%d] IceModel_PISM::run_timestep: ivals2.size()=%ld does not match pism_ivars.size()=%ld", pism_rank, ivals2.size(), pism_ivars.size());

    // Check Petsc types
    if (sizeof(double) != sizeof(PetscScalar)) {
        (*icebin_error)(-1, "PetscScalar must be same as double\n");
    }

    // Create a sorted order of our indices
    // (so we can eliminate duplicates via consolidate_by_perm)
    std::vector<int> perm(spsparse::sorted_perm(indices));

    // --------- Create g2_ix: PISM index of each value we received.
    // nvals: the LARGEST our resulting vector might have
    //     to be, if we don't have any duplicates.  In 
    // nconsolidated: the ACTUAL size of our resulting vector.

    // Remove duplicates (and sort too)
    blitz::Array<int,1> g2_ix(nvals);
    int nconsolidated = spsparse::consolidate_by_perm(indices, perm,
        indices, g2_ix, spsparse::DuplicatePolicy::REPLACE);

#if 0
    // No need for this, now that IceBin can do PISM-native indexing.
    // Convert to PISM indexing
    for (int ix=0; ix<nconsolidated; ++ix) {
        g2_ix(ix) = icebin_to_pism1d(g2_ix(ix));
    }
#endif

    // -----------------
    // Use the uniq-ified indices to fill pism_ivars[i] <-- ivals2[i]
    // pism_ivars are distributed (global) vectors.
    blitz::Array<PetscScalar,1> g2_y(nconsolidated);
    for (unsigned int i=0; i<contract[IceModel::INPUT].index.size(); ++i) {
        VarMeta &cf = contract[IceModel::INPUT].data[i];

        // Get matching input (val) and output (pism_var) variables
        IceModelVec2S *pism_var = pism_ivars[i];

        // Inputs specified in the contract are not (necessarily) attached
        // to any PISM var.  If they are not, just drop them on the ground.
        if (!pism_var) continue;

        // Consolidate the value array
        // Now we have a sparse vector, represented as (g2_ix, g2_y)
        g2_y = 0.;
        spsparse::consolidate_by_perm(indices, perm,
            ivals2[i], g2_y, spsparse::DuplicatePolicy::ADD);
//g2_y = 250.;

        // Put into a natural-ordering global distributed Petsc Vec
        ierr = VecSet(g2natural, cf.default_value); PISM_CHK(ierr, "run_timestep");
        ierr = VecSetValues(g2natural, nconsolidated, &g2_ix(0), &g2_y(0), INSERT_VALUES); PISM_CHK(ierr, "run_timestep");

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

printf("IceModel_PISM::run_timestep(): timestep_s = %g\n", timestep_s);

    // -------- Determine Dirichlet B.C. for ice sheet
    // This is done by taking the changes in the "borrowed" enthalpies
    // from the GCM, and applying them to the corresponding top layer in
    // the ice model.  The result is placed in surface->surface_temp.
    pism::icebin::IBSurfaceModel * const surface = ice_model->ib_surface_model();
    ice_model->construct_surface_temp(
        surface->icebin_deltah,
        contract[INPUT].at("deltah").default_value,
        timestep_s,
        surface->surface_temp);

    pism_in_nc->write(time_s);


    // =========== Run PISM for one coupling timestep
    // Time of last time we coupled
    printf("BEGIN ice_model->run_to(%f -> %f) %p\n",
        pism_grid->ctx()->time()->current(), time_s, ice_model.get());
    // See icebin::gpism::PISMIceModel::run_to()
    ice_model->run_to(time_s);
    printf("END ice_model->run_to()\n");


    if ((ice_model->mass_t() != time_s) || (ice_model->enthalpy_t() != time_s)) {
        (*icebin_error)(-1,
            "ERROR: PISM time (mass=%f, enthalpy=%f) doesn't match ICEBIN time %f", ice_model->mass_t(), ice_model->enthalpy_t(), time_s);
    }

    // ice_model->enthalpy_t() == time_s here
    ice_model->prepare_outputs(old_pism_time);
printf("pism_out_nc->write() after regular timestep\n");
    pism_out_nc->write(time_s);
//    ice_model->pism_out_state_nc->write(time_s);

    get_state(0);
    ice_model->reset_rate();

printf("Current time is pism: %f --> %f, ICEBIN: %f\n", old_pism_time, pism_grid->ctx()->time()->current(), time_s);

    printf("END IceModel_PISM::run_timestep(%f)\n", time_s);
}


/** Copies PISM->Icebin output variables from PISM variables to
the Icebin-supplied variables (on the root node).
@param mask Only do it for variables where (flags & mask) == mask.  Set to 0 for "all." */
void IceModel_PISM::get_state(unsigned int mask)
{
    printf("BEGIN IceModel_PISM::get_state: %ld\n", pism_ovars.size());
    VarSet const &ocontract(contract[IceModel::OUTPUT]);

    // Copy the outputs to the blitz arrays
    if (am_i_root()) allocate_ice_ovals_I();        // ROOT in PISM communicator
    for (unsigned int i=0; i<pism_ovars.size(); ++i) {
        if (!pism_ovars[i]) (*icebin_error)(-1,
            "IceModel_PISM: Contract output %s (modele_pism.cpp) is not linked up to a pism_ovar (MassEnergyBudget.cpp)", ocontract.index[i].c_str());

        VarMeta const &cf(ocontract.data[i]);
        if ((cf.flags & mask) != mask) continue;

        printf("IceModel_PISM::get_state(mask=%d) copying field %s\n", mask, cf.name.c_str());

        if (am_i_root()) {      // ROOT in PISM communicator

            // Check number of variables matches for output
            if (ice_ovals_I.size() != pism_ovars.size()) (*icebin_error)(-1,
                "[%d] IceModel_PISM::run_timestep: ice_ovals_I.size()=%ld does not match pism_ovars.size()=%ld", pism_rank, ice_ovals_I.size(), pism_ovars.size());

            // Reshape 1D blitz variable to 2D for use with PISM
            blitz::Array<double,2> oval2_xy(
                ibmisc::reshape<double,1,2>(ice_ovals_I[i], blitz::shape(ny(), nx())));
        
    
            // Get matching input (val) and output (pism_var) variables
            iceModelVec2S_to_blitz_xy(*pism_ovars[i], oval2_xy);    // Allocates oval2_xy if needed

        } else {
            blitz::Array<double,2> oval2_xy;    // dummy
            iceModelVec2S_to_blitz_xy(*pism_ovars[i], oval2_xy);
        }

        // Now send those data from the PISM root to the GCM root (MPI nodes)
        // (DUMMY for now, just make sure PISM and GCM have the same root)
        if (pism_root != coupler->gcm_params.gcm_root) (*icebin_error)(-1,
            "PISM and the GCM must share the same root!");
    }
    printf("END IceModel_PISM::get_state\n");
}


void IceModel_PISM::get_initial_state(double time_s)
{
    PetscErrorCode ierr;

//  double time_s = coupler->gcm_params.time_start_s;

    // Only prepare PISMIceModel outputs for things we need at init time.
    ice_model->prepare_initial_outputs(); PISM_CHK(ierr, "run_timestep");
printf("pism_out_nc->write() after get_initial_state\n");
    pism_out_nc->write(time_s); PISM_CHK(ierr, "run_timestep");
//    ice_model->pism_out_state_nc->write(time_s); PISM_CHK(ierr, "run_timestep");

    // Copy outputs to Icebin-supplied variables, only for INITIAL variables
printf("[%d] Calling get_state(%d)\n", pism_rank, contracts::INITIAL);
    get_state(contracts::INITIAL); PISM_CHK(ierr, "run_timestep");
}


}}
