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
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <icebin/GCMCoupler.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/contracts/contracts.hpp>
#include <spsparse/multiply_sparse.hpp>
#include <spsparse/sort.hpp>

#ifdef USE_PISM
#include <icebin/pism/IceModel_PISM.hpp>
#endif

using namespace spsparse;
using namespace ibmisc;
using namespace netCDF;

namespace icebin {


std::unique_ptr<IceCoupler> new_ice_coupler(NcIO &ncio, std::string vname,
    GCMCoupler const *_gcm_coupler, IceRegridder *_regridder)
{
    std::string vn(vname + ".info");
    auto info_v = get_or_add_var(ncio, vn, "int64", {});

    IceCoupler::Type type;
    get_or_put_att_enum(info_v, ncio.rw, "ice_coupler", type);

    std::unique_ptr<IceCoupler> ice_coupler;
    switch(type.index()) {
#if 0
        case IceCoupler::Type::DISMAL :
            ice_coupler.reset(new IceCoupler_DISMAL);
        break;
#endif
#ifdef USE_PISM
        case IceCoupler::Type::PISM :
            ice_coupler.reset(new gpism::IceCoupler_PISM);
        break;
#endif
        default :
            (*icebin_error)(-1,
                "Unknown IceCoupler::Type %s", type.str());
    }


    // Do basic initialization...
    ice_coupler->coupler = _coupler;
    ice_coupler->sheet = _sheet;
//    ice_coupler->ice_constants.init(&_coupler->ut_system);

    ice_coupler->gcm_per_ice_sheet_params =
        gcm_coupler->read_gcm_per_ice_sheet_params(ncio, vname);
    ice_coupler->ncread(ncio, vname);

    return ice_coupler;
}


IceCoupler::~IceCoupler() {}
// ==========================================================
static double const nan = std::numeric_limits<double>::quiet_NaN();

// ------------------------------------------------------------

/** Returns just the matrix part of regridding.  Also sets scale and
    correctA.  Local convenience function. */
inline SparseMatrix regrid_M(
IceRegridder &regridder,
std::string const &spec_name)
    { return std::move(regridder.regrid(spec_name, true, true)->M); }



// ==============================================================
void IceCoupler::set_start_time(
        ibmisc::time::tm const &time_base,
        double time_start_s)
{
    // Set up writers
    for (int io=0; io<2; ++io) {
        ice_coupler->writer[io].reset(new IceWriter(
            *this, contract[io],
            name() + (io == 0 ? "_in.nc" : "_out.nc")));
    }

    // This function is the last phase of initialization.
    // Only now can we assume that the contracts are fully set up.
#if 1
    // Print out the contract and var transformations
    std::cout << "========= Contract for " << name() << std::endl;
    std::cout << "---- GCM->Ice     Output Variables:" << std::endl;
    std::cout << contract[IceCoupler::INPUT];
    std::cout << "TRANSFORMATIONS:" << std::endl;
    std::cout << var_transformer[IceCoupler::INPUT];
    std::cout << "---- Ice->GCM     Output Variables:" << std::endl;
    std::cout << contract[IceCoupler::OUTPUT];
    std::cout << "TRANSFORMATIONS:" << std::endl;
    std::cout << var_transformer[IceCoupler::OUTPUT];
#endif


}

// ==============================================================

static std::array<std::string, 2> _writer_ofname = {"ice_model_in.nc", "ice_model_out.nc"};


/** 
@param do_run True if we are to actually run (otherwise just return ice_ovalsI from current state)
@return ice_ovalsI */
void IceCoupler::couple(
double time_s,
// Values from GCM, passed GCM -> Ice
ArraySparseParallelVectorsE const &gcm_ovalsE,
GCMCoupleOutput &out,    // Accumulate matrices here...
bool do_run)
{
    // Store regridding matrices for the last timestep, which we will
    // need to create ice_ivals
//    std::unique_ptr<WeightedSparse> IvEd0(std::move(IvEd));    // Dense A and E
//    std::unique_ptr<WeightedSparse> IvAd0(std::move(IvAd));

    // ========== Get Ice Inputs
    blitz::Array<double,2> ice_ivalsI({ndata(), contract[INPUT].size()});
    blitz::Array<double,2> ice_ovalsI({ndata(), contract[OUTPUT].size()});
    ice_ivalsI = 0;
    ice_ovalsI = 0;

    // Save old IvE for now
    auto IvE0(std::move(IvE));

    // ------------- Form ice_ivalsI
    // Densify gcm_ovals and IvE
    SparseSetT dimE;
    SparseMatrix IvEd0;
    dimE.add_sorted(IvE0.dim_begin(1), IvE0.dim_end(1));
    densify_one_dim(IvEd0, IvE0, dimE, 1);

    // Densify gcm_ovalsE --> gcm_ovalsEd
    // This should ONLY involve iE already mentioned in IvE0;
    // if not, there will be an exception.
    blitz::Array<double,2> gcm_ovalsEd(dimE.sparse_extent(), gcm_coupler->gcm_outputsE);
    gcm_ovalsEd = 0;
    for (size_t i=0; i<gcm_ovalsE.index.size(); ++i) {
        auto iE(gcm_ovalsE.index[i]);
        iEd = dimE.to_dense(iE);
        for (size_t ivar=0; ivar<gcm_ovalsE.values.size(); ++i) {
            gcm_ovals(iEd, ivar) += gcm_ovalsE.values[ivar](i);
        }
    }

    // Get the CSR sparse matrix to convert GCM outputs to ice model inputs
    auto scalars({
        std::make_pair("by_dt", 1.0 / ((itime - api->itime_last) * api->dtsrc)),
        std::make_pair("unit", 1.0)});
    CSRAndUnits icei_v_gcmo(var_transformer[INPUT].apply_scalars(scalars));

    // ice_ivalsEd_{jk} = icei_v_gcmo_{kl} * gcm_ovalsEd_{jl}
    // ice_ivalsI_{ik} = IvEd_{ij} * ice_ivalsEd_{jk}
    //       or:
    // ice_ivalsI_{ik} = IvEd_{ij} * icei_v_gcmo_{kl} * gcm_ovalsEd_{jl}
    //       where:
    // |i| = # ice grid cells (|I|)
    // |j| = # dense elevation grid cells (|Ed|)
    // |k| = # variables in ice_input
    // |l| = # variables in gcm_output
    //
    // (NOTE storage order; indices are row-major)

    // Regrid & combine to form ice_ivalsI
    for (auto iIvEd = IvEd.begin(); iIvEd != IvEd.end(); ++iIvEd) {
        auto ii(iIvEd.index(0));
        auto jj(iIvEd.index(1));
        auto IvEd_ij(iIvEd.val());

        // Transform units on the input while multiplying by M
        for (size_t kk = 0; kk < contracts[INPUT].size(); ++kk) {
            double ice_ivalsEd_jk = 0;
            std::vector<std::pair<int, double>> const &row(icei_v_gcmo.mat[kk]);
            for (auto rowk_iter=row.begin(); rowk_iter != row.end(); ++rowk_iter) {
                auto ll(rowk_iter->first);
                auto icei_v_gcmo_kl(row_iter->second);

                ice_ivalsEd_jk += icei_v_gcmo_kl * gcm_ovalsEd(jj, ll);
            }
            ice_ivalsI(ii, kk) += IvEd_ij * ice_ivalsEd_jk;
        }
    }

    // ========= Step the ice model forward
    if (writer[INPUT].get()) writer[INPUT]->write(time_s, ice_ivalsI);
    run_timestep(time_s, ice_ivalsI, ice_ovalsI, do_run);
    if (writer[OUTPUT].get()) writer[OUTPUT]->write(time_s, ice_ovalsI);

    // ========== Update regridding matrices
    auto elevI(get_elevI());    // blitz
    RegridMatrices rm(regridder);

    // Compute IvE (for next timestep)
    IvE = regrid_M(re, "IvE");
    dimE.add_sorted(IvE->M.dim_begin(1), IvE->M.dim_end(1));
    densify_one_dim(IvEd, IvE->M, dimE, 1);

    // Compute regrid matrices we need now
    auto EvI(regrid(re, "EvI"));
    auto AvI(regrid(re, "AvI"));
    auto AvE(regrid(re, "AvE"));
    copy(out.AvE1, AvE.M, AvE.dims[0], AvE.dims[1]);
    copy(out.wAvE1, AvE.weight, AvE.dims[0]);    // Blitz array (or should it be an Eigen column vector?)

    // Accumulate global things for all ice sheets....
    copy(out.E1vE0, EvI.M*IvE0, EvI.dims[0], dimE0);   // copy eigen matrix
    copy(out.elevE1, EvI * to_col_vector(elevI), EvI.dims[0]);   // copy eigen column vector

    // ========= Compute gcm_ivalsE = EvI * vt * ice_ovals
    CSRAndUnits gcmi_v_iceo(var_transformer[OUTPUT].apply_scalars(scalars));

    std::array<WeightedSparse *, GCMCoupler::GCMI::COUNT> XvIs;
    XvIs[GCMCoupler::GCMI::E] = &EvI;
    XvIs[GCMCoupler::GCMI::A] = &AvI;
    std::vector<double> vals(contract.size());

    // Do it once for _E variables and once for _A variables.
    for (int iEA=0; iEA < GCMCoupler::GCMI::COUNT; ++iEA) {
        VarSet &contract(gcm_coupler->gcm_inputs[iEA]);
        WeightedSparse &XvI_ws(*XvIs[iEA]);
        SparseMatrix XvI(to_spsparse(XvI_ws.M, XvI_ws.dims[0], XvI_ws.dims[1]));
        SparseParallelVectors &gcm_ivalsX(coupler_ret.gcm_ivals[iEA]);

        // gcm_ivalsX_{jn} = XvI_{ji} * gcmi_v_iceo_{nm} * ice_ovalsI_{im}
        //       where:
        // |i| = # ice grid cells (|I|)
        // |j| = # elevation/atmosphere grid cells (|X|)
        // |m| = # variables in ice_output
        // |n| = # variables in gcm_input

        // Do the multiplication
        for (auto iXvI(XvI->begin()); iXvI != XvI->end(); ++iXvI) {
            auto jj(iXvI.index(0));
            auto ii(iXvI.index(1));
            auto XvI_ji(iXvI.val());

            // Transform units on the input while multiplying by M
            for (size_t nn = 0; nn < contract.size(); ++nn) {
                double zval = 0;
                double gcm_ivalsI_in = 0;
                std::vector<std::pair<int, double>> const &row(gcmi_v_iceo.mat[nn]);
                for (auto rown_iter=row.begin(); rown_iter != row.end(); ++rown_iter) {
                    auto mm(rown_iter->first);
                    auto gcmi_v_iceo_nm(rown_iter->second);
                    gcm_ivalsX_in += gcmi_v_iceo_nm * ice_ovalsI(ii, mm);
                }
                vals[nn] = XvI_ji * gcm_ivalsX_in;
            }
            gcm_ivals.add(jj, vals);
        }
    }

    dimE0 = std::move(dimE);
}
