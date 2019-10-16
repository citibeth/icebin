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
    for (int nv : nvar) {
        gcm_ivalss_s.push_back(VectorMultivec(nv));
        gcm_ivalss_weight_s.push_back(std::vector<double>());
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

    ice_couplers.clear();
    for (size_t i=0; i < gcm_regridder->ice_regridders().size(); ++i) {
        std::string const &sheet_name(gcm_regridder->ice_regridders()[i]->name());

        // Create an IceCoupler corresponding to this IceSheet.
        std::unique_ptr<IceCoupler> ice_coupler(new_ice_coupler(ncio_config, vname, sheet_name, this));

        ice_couplers.push_back(std::move(ice_coupler));
    }


    printf("END GCMCoupler::ncread(%s)\n", grid_fname.c_str()); fflush(stdout);
}



void GCMCoupler::cold_start(
    ibmisc::Datetime _time_base,
    double _time_start_s)
{
    printf("BEGIN GCMCoupler::cold_start() nsheets=%ld\n", ice_couplers.size());
    time_base = _time_base;
    time_start_s = _time_start_s;
    time_unit = TimeUnit(&cal365, time_base, TimeUnit::SECOND);
    last_time_s = time_start_s;

    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        auto &ice_coupler(ice_couplers[sheetix]);

        // NOTE: Contract sizes are only needed at runtime, not allocate() time.
//        contracts::setup(*this, *ice_coupler);    // where does this go w.r.t ncread() and upate?

        // Dynamic ice model is instantiated here...
        ice_coupler->cold_start(time_base, time_start_s);

        ice_coupler->print_contracts();
    }

    printf("END GCMCoupler::cold_start()\n");
}
// ------------------------------------------------------------

static void ncwrite_dense_VectorMultivec(
    NcGroup * const nc,
    VectorMultivec const *vecs,
    VarSet const *contract,
    ibmisc::Indexing const *indexing,
    std::string const &vname_base)
{
printf("BEGIN ncwrite_dense_VectorMultivec()\n");

    // im,jm,ihc  0-based
    long nE = indexing->extent();
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

    // Go through each variable...
    unsigned int nvar = vecs->nvar;
    for (unsigned int ivar=0; ivar < nvar; ++ivar) {
        // Fill our dense var
        denseE = nan;
        for (unsigned int i=0; i<vecs->index.size(); ++i) {
            auto iE(vecs->index[i]);
            if (iE >= nE) (*icebin_error)(-1,
                "Index out of range: %ld vs. %ld", (long)iE, (long)nE);
            denseE(iE) = vecs->vals[i*nvar + ivar];
        }

        // Store in the netCDF variable
        NcVar ncvar(nc->getVar(vname_base + (*contract)[ivar].name));
        ncvar.putVar(startp, countp, denseE.data());
    }

printf("END ncwrite_dense()\n");
}

/** Densifies the sparse vectors, and then writes them out to netCDF */
static void ncio_dense(
    NcIO &ncio,
    VectorMultivec const &vecs,    
    VarSet const &contract,
    ibmisc::Indexing const &indexing,
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
        ncio.nc, &vecs, &contract, &indexing, vname_base));

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
        ncio.nc, &E_s, &indexing, vname));

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

printf("indexing extent1 = %ld\n", indexing->extent());

        ncio_dense(ncio, out.gcm_ivalss_s[iAE], gcm_inputs[iAE],
            *indexing, vname_base);

        auto dims(get_or_add_dims(ncio,
            {"gcm_ivalss_weight_"+IndexAE_labels[iAE]+"_len"},
            {out.gcm_ivalss_weight_s[iAE].size()}));
        ncio_vector(ncio, out.gcm_ivalss_weight_s[iAE], false,
            "gcm_ivalss_weight_"+IndexAE_labels[iAE]+"_s", "double", dims);
    }

//    ncio_spsparse(ncio, out.E1vE0_unscaled, false, vname_base+"E1vE0_unscaled");
    out.E1vE0_unscaled.ncio(ncio, vname_base+"E1vE0_unscaled");
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
        gcm_regridder->indexingE, vname_base);
printf("END GCMCoupler::ncio_gcm_output(%s)\n", vname_base.c_str());
}
// ------------------------------------------------------------
// ======================================================================

}       // namespace
