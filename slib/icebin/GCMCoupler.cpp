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

#include <mpi.h>        // Intel MPI wants to be first
#include <functional>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/ncfile.hpp>
#include <ibmisc/memory.hpp>
#include <ibmisc/blitz.hpp>
#include <ibmisc/datetime.hpp>
#include <icebin/GCMCoupler.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/contracts/contracts.hpp>
#include <icebin/e1ve0.hpp>
#include <spsparse/netcdf.hpp>

#ifdef USE_PISM
#include <icebin/pism/IceCoupler_PISM.hpp>
#endif

using namespace spsparse;
using namespace ibmisc;
using namespace netCDF;

namespace icebin {

static double const nan = std::numeric_limits<double>::quiet_NaN();

// ==========================================================


GCMParams::GCMParams(MPI_Comm _gcm_comm, int _gcm_root) :
    gcm_comm(_gcm_comm),
    gcm_root(_gcm_root),
    world(gcm_comm, boost::mpi::comm_attach)
{
    MPI_Comm_rank(gcm_comm, &gcm_rank);
}

// ==========================================================
/** @param nvar Array specifying number of variables for each segment. */
GCMInput::GCMInput(std::vector<int> const &nvar) {
    for (size_t ix=0; ix<nvar.size(); ++ix) {
        gcm_ivalss_s.push_back(VectorMultivec(nvar[ix]));
    }
}

std::vector<int> GCMInput::nvar() const
{
    std::vector<int> ret;
    ret.reserve(gcm_ivalss_s.size());
    for (size_t i=0; i<gcm_ivalss_s.size(); ++i) ret.push_back(gcm_ivalss_s[i].nvar);
    return ret;
}


// ==========================================================

GCMCoupler::GCMCoupler(Type _type, GCMParams &&_gcm_params) :
    type(_type),
    gcm_params(std::move(_gcm_params)),
    ut_system(""),
    use_smb(true)
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


/** Produces a date string in format YYMMDD */
std::string GCMCoupler::sdate(double time_s) const
{
    ibmisc::Datetime dt(time_unit.to_datetime(time_s));
    return (boost::format
        ("%04d%02d%02d") % dt[0] % dt[1] % dt[2]).str();
}


/** Static-cast a unique_ptr while moving it to another unique_ptr */
template<typename D, typename B>
void static_move(std::shared_ptr<D> &dest, std::unique_ptr<B>& base)
{
    dest.reset(static_cast<D*>(base.release()));
}


/** @param nc The IceBin configuration file */
void GCMCoupler::_ncread(
    ibmisc::NcIO &ncio_config,
    std::string const &vname)        // comes from this->gcm_params
{
    auto config_info(get_or_add_var(ncio_config, vname + ".info", "int", {}));
    std::string grid_fname;
    get_or_put_att(config_info, ncio_config.rw, "grid", grid_fname);
    get_or_put_att(config_info, ncio_config.rw, "output_dir", output_dir);
    get_or_put_att(config_info, ncio_config.rw, "use_smb", &use_smb, 1);

    printf("BEGIN GCMCoupler::ncread(%s)\n", grid_fname.c_str()); fflush(stdout);

//    bool rw_full = am_i_root();

    // Load the MatrixMaker (filtering by our domain, of course)
    // Also load the ice sheets
    {
        std::unique_ptr<GCMRegridder_Standard> gcmr(new GCMRegridder_Standard());
        NcIO ncio_grid(grid_fname, NcFile::read);
        gcmr->ncio(ncio_grid, vname);
        static_move(gcm_regridder, gcmr);    // Move gcm_regridder <- gcm
    }

    std::cout << "========= GCM Constants" << std::endl;
    std::cout << gcm_constants;
    std::cout << "========= GCM Outputs" << std::endl;
    std::cout << gcm_outputsE;
    std::cout << "========= GCM Inputs" << std::endl;
    for (size_t i=0; i<gcm_inputs.size(); ++i) {
        std::cout << "--------------- " << i << " " << gcm_inputs_grid[i] << std::endl;
        std::cout << gcm_inputs_grid[i] << ':' << gcm_inputs[i];
    }

    // Compute static info about combined exchange grids
    size_t const nsheet = gcm_regridder->ice_regridders().size();
    ice_couplers.clear();
    basesX = std::vector<long>{0};  // Base of first exchange grid is 0
    areaX.clear();
    for (auto &ice_regridder : gcm_regridder->ice_regridders()) {
        // Create an IceCoupler corresponding to this IceSheet.
        ice_couplers.push_back(new_ice_coupler(
            ncio_config, vname, ice_regridder->name(), this));

        // Add to basesX and areaX (X = combined exchange grid for all ice sheets)
        auto const &aexgrid(ice_regridder->aexgrid);    // Exchange grid
        basesX.push_back(basesX.back() + aexgrid.sparse_extent());
        for (int iXd=0; iXd<aexgrid.dense_extent(); ++iXd) {
            areaX.push_back(aexgrid.native_area(iXd));
        }
    }

    printf("END GCMCoupler::ncread(%s)\n", grid_fname.c_str()); fflush(stdout);
}



void GCMCoupler::model_start(
    bool cold_start,
    ibmisc::Datetime _time_base,
    double _time_start_s)
{
    printf("BEGIN GCMCoupler::model_start() nsheets=%ld\n", ice_couplers.size());
    time_base = _time_base;
    time_start_s = _time_start_s;
    time_unit = TimeUnit(&cal365, time_base, TimeUnit::SECOND);
    timespan = std::array<double,2>{-1,time_start_s};

    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        auto &ice_coupler(ice_couplers[sheetix]);

        // NOTE: Contract sizes are only needed at runtime, not allocate() time.
//        contracts::setup(*this, *ice_coupler);    // where does this go w.r.t ncread() and upate?

        // Dynamic ice model is instantiated here...
        ice_coupler->model_start(cold_start, time_base, time_start_s);
        ice_coupler->print_contracts();
    }

    printf("END GCMCoupler::model_start()\n");
}
// ------------------------------------------------------------

static void ncwrite_dense_VectorMultivec(
    NcGroup * const nc,
    VectorMultivec const *vecs,
    VarSet const *contract,
    ibmisc::Indexing const *indexing,
    ibmisc::UTSystem const *ut_system,
    std::string const &vname_base)
{
    // NOTE: This code is in E (elevation grid); but it works just as
    // well for A (atmosphere grid)
printf("BEGIN ncwrite_dense_VectorMultivec(nnz=%ld)\n", vecs->index.size());

    // im,jm,ihc  0-based
    long nE = indexing->extent();
    blitz::Array<double,1> scaleE(nE);
    blitz::Array<double,1> denseE(nE);    // All dims

    // NetCDF needs dimensions in stride-descending order
    std::vector<size_t> startp;
    std::vector<size_t> countp;
    startp.push_back(0);    // Time dimension
    countp.push_back(1);
    for (int dimi : indexing->indices()) {
        startp.push_back(0);
        countp.push_back((*indexing)[dimi].extent);
    }

    vecs->to_dense_scale(scaleE);

    // Go through each variable...
    unsigned int nvar = vecs->nvar;
    for (unsigned int ivar=0; ivar < nvar; ++ivar) {
        vecs->to_dense(ivar, scaleE, nan, denseE);

        // Convert to NC units for output file
        denseE *= (*contract)[ivar].nc_factor(*ut_system);

        // Store in the netCDF variable
        NcVar ncvar(nc->getVar(vname_base + (*contract)[ivar].name));
        ncvar.putVar(startp, countp, denseE.data());
    }

printf("END ncwrite_dense_VectorMultivec()\n");
}

/** Densifies the sparse vectors, and then writes them out to netCDF */
static void ncio_dense(
    NcIO &ncio,
    VectorMultivec const &vecs,    
    VarSet const &contract,
    ibmisc::Indexing const &indexing,
    UTSystem const &ut_system,
    std::string const &vname_base = "")
{
printf("BEGIN ncio_dense1('%s')\n", ncio.fname.c_str());
    if (ncio.rw != 'w') (*icebin_error)(-1,
        "ncio_dense(VectorMultivec) only writes, no read.");

    if (indexing.extent() == 0) (*icebin_error)(-1,
        "indexing.extent() == 0");

    // Create the variables
    NcDimSpec dim_spec({"time"}, {1});
    append(dim_spec, indexing);
    contract.ncdefine(ncio, dim_spec.to_dims(ncio), vname_base);

    ncio.add("ncio_dense::" + vname_base, std::bind(ncwrite_dense_VectorMultivec,
        ncio.nc, &vecs, &contract, &indexing, &ut_system, vname_base));

printf("END ncio_dense()\n");
}
// --------------------------------------------------------------
#if 0	// Not needed
static void ncwrite_dense_TupleList1(
    NcGroup * const nc,
    TupleListLT<1> const *E_s,
    ibmisc::Indexing const *indexing,
    std::string const &vname)
{
printf("BEGIN ncwrite_dense()\n");

    // --------- Convert sparse to dense (1-D indexing)
    long nE = indexing->extent();
    blitz::Array<double,1> denseE(nE);    // All dims
    for (auto ii = E_s->begin(); ii != E_s->end(); ++ii) {
        // 1-D indexing here
        auto iE(ii->index(0));
        denseE(iE) = ii->value();
    }

    // NetCDF needs dimensions in stride-descending order
    std::vector<size_t> startp;
    std::vector<size_t> countp;
    startp.push_back(0);    // Time dimension
    countp.push_back(1);
    for (int dimi : indexing->indices()) {
        startp.push_back(0);
        countp.push_back((*indexing)[dimi].extent);
    }

    // Store in the netCDF variable
    NcVar ncvar(nc->getVar(vname));
    ncvar.putVar(startp, countp, denseE.data());

printf("END ncwrite_dense()\n");
}

/** Densifies the sparse vectors, and then writes them out to netCDF */
static void ncio_dense(
    NcIO &ncio,
    TupleListLT<1> const &E_s,
    ibmisc::Indexing const &indexing,
    ibmisc::UTSystem const &ut_system,
    std::string const &vname)
{
printf("BEGIN ncio_dense2('%s')\n", ncio.fname.c_str());
    if (ncio.rw != 'w') (*icebin_error)(-1,
        "ncio_dense(TupleListT<1>) only writes, no read.");

    // Create the variables
    NcDimSpec dim_spec({"time"}, {1});
    append(dim_spec, indexing);
    get_or_add_var(ncio, vname, "double", dim_spec.to_dims(ncio));

    ncio.add("ncio_dense::" + vname, std::bind(ncwrite_dense_TupleList1,
        ncio.nc, &E_s, &indexing, &ut_system, vname));

printf("END ncio_dense()\n");
}
#endif
// ------------------------------------------------------------
/** Top-level ncio() to log output from coupler. (coupler->GCM) */
void GCMCoupler::ncio_gcm_input(NcIO &ncio,
    GCMInput &out,        // All MPI ranks included here
    std::array<double,2> &timespan,    // timespan[1] = current time_s
    ibmisc::TimeUnit const &time_unit,
    std::string const &vname_base)
{
    ibmisc::ncio_timespan(ncio, timespan, time_unit, vname_base + "timespan");

    for (int iAE=0; iAE<(int)IndexAE::COUNT; ++iAE) {
        ibmisc::Indexing const *indexing;
        switch(iAE) {
            case (int)IndexAE::A:
            case (int)IndexAE::ATOPO:
                indexing = &gcm_regridder->indexing(GridAE::A);
            break;
            case (int)IndexAE::E:
            case (int)IndexAE::ETOPO:
                indexing = &gcm_regridder->indexing(GridAE::E);
            break;
            default:
                (*icebin_error)(-1, "Unknown iAE=%d; this switch statement needs to be expanded?", iAE);
            break;
        }

        // The indexing is used from gcm_regridder.  This will result in nhc_ice
        ncio_dense(ncio, out.gcm_ivalss_s[iAE], gcm_inputs[iAE],
            *indexing, ut_system, vname_base);

    }

//    ncio_spsparse(ncio, out.E1vE0_unscaled, false, vname_base+"E1vE0_unscaled");
    out.E1vE0c.ncio(ncio, vname_base+"E1vE0c");
}

/** Top-level ncio() to log input to coupler. (GCM->coupler) */
void GCMCoupler::ncio_gcm_output(NcIO &ncio,
    VectorMultivec const &gcm_ovalsE,
    std::array<double,2> &timespan,    // timespan[1] = current time_s
    ibmisc::TimeUnit const &time_unit,
    std::string const &vname_base)
{
printf("BEGIN GCMCoupler::ncio_gcm_output('%s' '%s')\n", ncio.fname.c_str(), vname_base.c_str());
    ncio_timespan(ncio, timespan, time_unit, vname_base + "timespan");
    ncio_dense(ncio, gcm_ovalsE, gcm_outputsE,
        gcm_regridder->indexingE, ut_system, vname_base);
printf("END GCMCoupler::ncio_gcm_output(%s)\n", vname_base.c_str());
}
// ------------------------------------------------------------
// ------------------------------------------------------------
// TODO: Add to linear/eigen.hpp
/** Converts a linear::Weighted_Eigen to sparse indexing
@param dims set to E.dims to sparsify all dimensions; or set one to nullptr if no sparsify needed there. */
std::unique_ptr<linear::Weighted_Eigen> sparsify(
    linear::Weighted_Eigen const &E,
    std::array<SparsifyTransform,2> const &transforms = {SparsifyTransform::TO_SPARSE, SparsifyTransform::TO_SPARSE},
//    std::array<SparseSetT *,2> const &dims,    // Determines what to sparsify; E.dims or nullptr in each slot
    std::array<int,2> extent = std::array<int,2>{-1,-1})            // Extent of resulting matrix (can be taken from E.dims if available)
{
    for (int i=0; i<2; ++i) if (extent[i] < 0) extent[i] = E.dims[i]->sparse_extent();
    std::unique_ptr<linear::Weighted_Eigen> S(
        new linear::Weighted_Eigen(E.dims, E.conservative));

    // Sparsify the matrix
    {
        TupleListT<2> SM_t;    // Matrix in sparse indexing, tuples
        spcopy(
            accum::sparsify(transforms, accum::in_index_type<int>(), E.dims,
            accum::ref(SM_t)),
            *E.M);

        S->M.reset(new EigenSparseMatrixT(
            extent[0], extent[1]));
        S->M->setFromTriplets(SM_t.begin(), SM_t.end());
    }

    // Sparsify the weights
    if (transforms[0] == SparsifyTransform::ID) {
        S->wM.reference(E.wM);
    } else {
        spcopy(
            accum::to_sparse(std::array<SparseSetT*,1>{E.dims[0]},
            accum::blitz_new(S->wM)),
            E.wM);
    }

    if (transforms[1] == SparsifyTransform::ID) {
        S->Mw.reference(E.Mw);
    } else {
        spcopy(
            accum::to_sparse(std::array<SparseSetT*,1>{E.dims[1]},
            accum::blitz_new(S->Mw)),
            E.Mw);
    }

    // Clear dims, which may not be saved...
    S->dims[0] = nullptr;
    S->dims[1] = nullptr;

    return S;
}
// -------------------------------------------------
GCMInput GCMCoupler::couple(
double time_s,        // Simulation time [s]
VectorMultivec const &gcm_ovalsE,
bool run_ice)    // if false, only initialize
{
    timespan = std::array<double,2>{timespan[1], time_s};

printf("BEGIN GCMCoupler::couple(time_s=%g, run_ice=%d)\n", time_s, run_ice);
    // ------------------------ Most MPI Nodes
    if (!gcm_params.am_i_root()) {
        GCMInput out({0,0,0,0});
        for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
            auto &ice_coupler(ice_couplers[sheetix]);
            ice_coupler->couple(timespan, gcm_ovalsE, out.gcm_ivalss_s, run_ice);
        }
        return out;
    }

    // ----------------------- Root MPI Node

    // -------- Figure out our calendar day to format filenames
    if (gcm_params.icebin_logging) {
        std::string fname = "gcm-out-" + this->sdate(time_s) + ".nc";
        NcIO ncio(fname, 'w');
        ncio_gcm_output(ncio, gcm_ovalsE, timespan,
            time_unit, "");
        ncio();
    }

    // ---------- Initialize output: A,E,Atopo,Etopo
    std::vector<int> nvars;
    for (auto &gcmi : gcm_inputs) nvars.push_back(gcmi.size());
    GCMInput out(nvars);

    // ---------- Run per-ice-sheet couplers
    {
        std::vector<SparseSetT const *> dimE1s;
        std::vector<std::unique_ptr<linear::Weighted_Eigen>> XuE1s;
        for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
            auto &ice_coupler(ice_couplers[sheetix]);    // IceCoupler
            IceRegridder const *ice_regridder = ice_coupler->ice_regridder;

            IceCoupler::CoupleOut iout(ice_coupler->couple(
                timespan, gcm_ovalsE, out.gcm_ivalss_s, run_ice));
            dimE1s.push_back(iout.dimE);
            XuE1s.push_back(sparsify(*iout.XuE,
                std::array<SparsifyTransform,2>{
                    SparsifyTransform::ID,
                    SparsifyTransform::TO_SPARSE},
                std::array<int,2>{ice_regridder->nX(), -1}));

        }

        // --------- Compute E1vE0
        if (run_ice) {
            out.E1vE0c = e1ve0::compute_E1vE0c(
                XuE1s, XuE0s,
                gcm_regridder->nE(), areaX);
        }
        XuE0s = std::move(XuE1s);    // save state between timesteps
    }

    return out;

}
// ------------------------------------------------------------
// ======================================================================

}       // namespace
