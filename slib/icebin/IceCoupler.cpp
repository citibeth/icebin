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
#include <type_traits>

#include <spsparse/blitz.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <ibmisc/array.hpp>
#include <ibmisc/string.hpp>    // string_printf()
#include <icebin/GCMCoupler.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/contracts/contracts.hpp>
#include <spsparse/eigen.hpp>
#include <spsparse/blitz.hpp>

#ifdef USE_PISM
#include <icebin/pism/IceCoupler_PISM.hpp>
#endif

using namespace spsparse;
using namespace ibmisc;
using namespace netCDF;
//using namespace std;


namespace icebin {

std::unique_ptr<IceCoupler> new_ice_coupler(NcIO &ncio,
    std::string const &vname, std::string const &sheet_name,
    GCMCoupler const *_gcm_coupler, IceRegridder *_ice_regridder)
{
    std::string vname_sheet(vname + "." + sheet_name);
    auto info_v = get_or_add_var(ncio, vname_sheet + ".info", "int64", {});

    IceCoupler::Type type;
    get_or_put_att_enum(info_v, ncio.rw, "ice_coupler", type);

    std::unique_ptr<IceCoupler> self;
    switch(type.index()) {
#if 0
        case IceCoupler::Type::DISMAL :
            self.reset(new IceCoupler_DISMAL);
        break;
#endif
#ifdef USE_PISM
        case IceCoupler::Type::PISM :
            self.reset(new gpism::IceCoupler_PISM);
        break;
#endif
        default :
            (*icebin_error)(-1,
                "Unknown IceCoupler::Type %s", type.str());
    }


    // Do basic initialization...
    self->_name = sheet_name;
    self->gcm_coupler = _gcm_coupler;
    self->ice_regridder = _ice_regridder;
    if (_ice_regridder) self->_nI = _ice_regridder->nI();

//    self->ice_constants.init(&_coupler->ut_system);

    self->ncread(ncio, vname_sheet);

    return self;
}


IceCoupler::~IceCoupler() {}
// ==========================================================
static double const nan = std::numeric_limits<double>::quiet_NaN();

// ==============================================================
void IceCoupler::cold_start(
        ibmisc::Datetime const &time_base,
        double time_start_s)
{
    // Set up writers
    if (gcm_coupler->am_i_root()) {
        for (int io=0; io<2; ++io) {    // INPUT / OUTPUT
            auto fname(
                boost::filesystem::path(output_dir) /
                (name() + (io == 0 ? "_in.nc" : "_out.nc")));
            this->writer[io].reset(new IceWriter(
                this, &contract[io], fname.string()));
        }
    }
}

/** Print summary info of the contracts to STDOUT. */
void IceCoupler::print_contracts()
{
    // This function is the last phase of initialization.
    // Only now can we assume that the contracts are fully set up.
#if 1
    // Print out the contract and var transformations
    std::cout << "========= Contract for " << name() << std::endl;
    std::cout << "---- Ice <- GCM     Output Variables:" << std::endl;
    std::cout << contract[IceCoupler::INPUT];
    std::cout << "TRANSFORMATIONS:" << std::endl;
    std::cout << var_trans_inE;
    std::cout << "---- GCM <- Ice     Output Variables:" << std::endl;
    std::cout << contract[IceCoupler::OUTPUT];
    std::cout << "TRANSFORMATIONS GCM-A <- Ice:" << std::endl;
    std::cout << var_trans_outAE[GridAE::A];
    std::cout << "TRANSFORMATIONS GCM-E <- Ice:" << std::endl;
    std::cout << var_trans_outAE[GridAE::E];
#endif


}

// ==============================================================

static std::array<std::string, 2> _writer_ofname = {"ice_model_in.nc", "ice_model_out.nc"};


static Eigen::VectorXd to_col_vector(blitz::Array<double,1> const &vec)
{
    Eigen::VectorXd ret(vec.extent(0));
    for (int i=vec.lbound(0), j=0; i <= vec.ubound(0); ++i,++j) ret(j) = vec(i);
    return ret;
}






/** 
@param do_run True if we are to actually run (otherwise just return ice_ovalsI from current state)
@return ice_ovalsI */
void IceCoupler::couple(
double time_s,
// Values from GCM, passed GCM -> Ice
VectorMultivec const &gcm_ovalsE_s,
GCMInput &out,    // Accumulate matrices here...
bool run_ice)
{
    if (!gcm_coupler->am_i_root()) {
        blitz::Array<double,1> elevI(get_elevI());    // blitz

        // Allocate dummy variables, even though they will only be set on root
        blitz::Array<double,2> ice_ivalsI(nI(), contract[INPUT].size());
        blitz::Array<double,2> ice_ovalsI(nI(), contract[OUTPUT].size());
        ice_ivalsI = 0;
        ice_ovalsI = 0;
        run_timestep(time_s, ice_ivalsI, ice_ovalsI, run_ice);
        return;
    }

    printf("BEGIN IceCoupler::couple(%s)\n", name().c_str());

    // This requires MPI...
    blitz::Array<double,1> elevI(get_elevI());    // blitz

    // ========== Get Ice Inputs
    blitz::Array<double,2> ice_ivalsI(nI(), contract[INPUT].size());
    blitz::Array<double,2> ice_ovalsI(nI(), contract[OUTPUT].size());
    ice_ivalsI = 0;
    ice_ovalsI = 0;

printf("ice_coupler BB1\n");
    // E_s = Elevation grid (sparse indices)
    // E0 = Elevation grid @ beginning of timestep (dense indices)
    // E1 = Elevation grid @ end of timestep (dense indices)
    // All matrices and vectors assumed w/ densified indexing
    // Except _s ending means they use sparse indexing.

    // ------------- Form gcm_ovalsE0
    // Densify gcm_ovalsE_s --> gcm_ovalsE
    // This should ONLY involve iE already mentioned in IvE0;
    // if not, there will be an exception.
    blitz::Array<double,2> gcm_ovalsE0(dimE0.dense_extent(), gcm_coupler->gcm_outputsE.size());
    gcm_ovalsE0 = 0;
    for (size_t i=0; i<gcm_ovalsE_s.index.size(); ++i) {
        auto iE_s(gcm_ovalsE_s.index[i]);
        int iE0 (dimE0.to_dense(iE_s));
        for (int ivar=0; ivar<gcm_ovalsE_s.size(); ++i) {
            gcm_ovalsE0(iE0, ivar) += gcm_ovalsE_s.val(ivar, i);
        }
    }

printf("ice_coupler BB2\n");
    // Get the CSR sparse matrix to convert GCM outputs to ice model inputs
    std::vector<std::pair<std::string, double>> scalars({
        std::make_pair("by_dt", 1.0 / (time_s - gcm_coupler->last_time_s)),
        std::make_pair("unit", 1.0)});
    CSRAndUnits icei_v_gcmo(var_trans_inE.apply_scalars(scalars));

    // ice_ivalsE0_{jk} = icei_v_gcmo_{kl} * gcm_ovalsE0_{jl}
    // ice_ivalsI_{ik} = IvE0_{ij} * ice_ivalsE0_{jk}
    //       or:
    // ice_ivalsI_{ik} = IvE0_{ij} * icei_v_gcmo_{kl} * gcm_ovalsE0_{jl}
    //       where:
    // |i| = # ice grid cells (|I|)
    // |j| = # dense elevation grid cells (|E0|)
    // |k| = # variables in ice_input
    // |l| = # variables in gcm_output
    //
    // (NOTE storage order; indices are row-major)

    // ------------- Form ice_ivalsI

printf("ice_coupler BB3\n");
    // Regrid & combine to form ice_ivalsI
    if (run_ice) {
    for (auto iIvE0(begin(*IvE0)); iIvE0 != end(*IvE0); ++iIvE0) {
        auto ii(iIvE0->index(0));
        auto jj(iIvE0->index(1));
        auto IvE0_ij(iIvE0->value());

        // Transform units on the input while multiplying by M
        for (int kk = 0; kk < contract[INPUT].size(); ++kk) {
            double ice_ivalsE_jk = 0;
            std::vector<std::pair<int, double>> const &row(icei_v_gcmo.mat[kk]);
            for (auto rowk_iter=row.begin(); rowk_iter != row.end(); ++rowk_iter) {
                int ll(rowk_iter->first);
                double icei_v_gcmo_kl(rowk_iter->second);
                ice_ivalsE_jk += icei_v_gcmo_kl * gcm_ovalsE0(jj, ll);
            }
            ice_ivalsI(ii, kk) += IvE0_ij * ice_ivalsE_jk;
        }
    }}

    // ========= Step the ice model forward
printf("ice_coupler BB4\n");
    if (writer[INPUT].get()) writer[INPUT]->write(time_s, ice_ivalsI);
printf("ice_coupler BB5\n");
    run_timestep(time_s, ice_ivalsI, ice_ovalsI, run_ice);
printf("ice_coupler BB6\n");
    if (writer[OUTPUT].get()) writer[OUTPUT]->write(time_s, ice_ovalsI);
printf("ice_coupler BB7\n"); fflush(stdout);

    // ========== Update regridding matrices
printf("ice_coupler BB7.1\n"); fflush(stdout);
    ice_regridder->set_elevI(elevI);
printf("ice_coupler BB7.2\n"); fflush(stdout);
    RegridMatrices rm(ice_regridder);
printf("ice_coupler BB7.3\n"); fflush(stdout);

    SparseSetT dimA1;
    SparseSetT dimE1;
    // A SparseSet that is identity for the entire range of I
    auto dimI(id_sparse_set<SparseSetT>(nI()));

printf("ice_coupler BB8\n"); fflush(stdout);
    // ---- Update AvE1 matrix and weights (global for all ice sheets)
    {
        auto AvE1(rm.regrid("AvE", {&dimA1, &dimE1}, true, true));

        spcopy(
            accum::to_sparse(AvE1->dims,
            accum::ref(out.AvE1_s)),
            *AvE1->M);

        TupleListLT<1> &out_wAvE1_s(out.wAvE1_s);    // Get around bug in gcc@4.9.3
        spcopy(
            accum::to_sparse(make_array(AvE1->dims[0]),
            accum::ref(out_wAvE1_s)),
            AvE1->weight);    // blitz::Array<double,1>
    }
printf("ice_coupler BB9.0\n");


    // ------ Update E1vE0 translation between old and new elevation classes
    //        (global for all ice sheets)
    auto E1vI(rm.regrid("EvI", {&dimE1, &dimI}, true, true));

    // Don't do this on the first round, since we don't yet have an IvE0
    if (run_ice) {
        TupleListLT<2> &out_E1vE0_s(out.E1vE0_s);
        EigenSparseMatrixT E1vE0(*E1vI->M * *IvE0);
        spcopy(
            accum::to_sparse(make_array(&dimE1, &dimE0),
            accum::ref(out_E1vE0_s)),
            E1vE0);
    }

printf("ice_coupler BB9.1\n"); fflush(stdout);
auto IvE1x(rm.regrid("IvE", {&dimI, &dimE1}, true, true));

printf("ice_coupler BB10\n"); fflush(stdout);


    // ------ Update orography (global for all ice sheets)
    TupleListLT<1> &out_elevE1_s(out.elevE1_s);
    Eigen::VectorXd elevE1(*E1vI->M * to_col_vector(elevI));
    spcopy(
        accum::permute(accum::in_rank<2>(), {0},   // copy eigen column vector
        accum::to_sparse(make_array(&dimE1),
        accum::ref(out_elevE1_s))),
        elevE1);

printf("ice_coupler BB11.0\n"); 
    // ========= Compute gcm_ivalsE = EvI * vt * ice_ovals

    auto A1vI(rm.regrid("AvI", {&dimA1, &dimI}, true, true));

    std::array<WeightedSparse * const, GridAE::count> AE1vIs {&*A1vI, &*E1vI};

printf("ice_coupler BB11.1\n");
auto IvE1x1(rm.regrid("IvE", {&dimI, &dimE1}, true, true));
printf("ice_coupler BB11.2\n"); 

// Somewhere in this #if block is corrupting memory
#if 1

    // Do it once for _E variables and once for _A variables.
    for (int iAE=0; iAE < GridAE::count; ++iAE) {
        CSRAndUnits gcmi_v_iceo(var_trans_outAE[iAE].apply_scalars(scalars));

        VarSet const &contract(gcm_coupler->gcm_inputsAE[iAE]);
        std::vector<double> vals(contract.size());
printf("ice_coupler BB11.3 [%d] vals.size() = %ld\n", iAE, vals.size());
        WeightedSparse &X1vI_ws(*AE1vIs[iAE]);
        VectorMultivec &gcm_ivalsX_s(out.gcm_ivalsAE[iAE]);

        // gcm_ivalsX_{jn} = X1vI_{ji} * gcmi_v_iceo_{nm} * ice_ovalsI_{im}
        //       where:
        // |i| = # ice grid cells (|I|)
        // |j| = # elevation/atmosphere grid cells (|X|)
        // |m| = # variables in ice_output
        // |n| = # variables in gcm_input

        // Do the multiplication
        for (auto iX1vI(begin(*X1vI_ws.M)); iX1vI != end(*X1vI_ws.M); ++iX1vI) {
            auto jj(iX1vI->index(0));
            auto ii(iX1vI->index(1));
            auto X1vI_ji(iX1vI->value());
//printf("X1vI[%d,%d] = %f\n",jj,ii,X1vI_ji);

            // Transform units on the input while multiplying by M
            for (size_t nn = 0; nn < contract.size(); ++nn) {
                double zval = 0;
                double gcm_ivalsX_in = 0;
                std::vector<std::pair<int, double>> const &row(gcmi_v_iceo.mat[nn]);
                for (auto rown_iter=row.begin(); rown_iter != row.end(); ++rown_iter) {
                    int const &mm(rown_iter->first);
                    double const &gcmi_v_iceo_nm(rown_iter->second);
                    gcm_ivalsX_in += gcmi_v_iceo_nm * ice_ovalsI((int)ii, (int)mm);
                }
                vals[nn] = X1vI_ji * gcm_ivalsX_in;
            }
            gcm_ivalsX_s.add(jj, vals);
        }
    }
#endif

printf("ice_coupler BB12\n");
    // Compute IvE (for next timestep)
    auto IvE1(rm.regrid("IvE", {&dimI, &dimE1}, true, true));
    dimE0 = std::move(dimE1);
    IvE0 = std::move(IvE1->M);
printf("ice_coupler BB13\n");

    printf("END IceCoupler::couple(%s)\n", name().c_str());
}
// =======================================================
/** Specialized init signature for IceWriter */
IceWriter::IceWriter(
    IceCoupler const *_ice_coupler,
    VarSet const *_contract,
    std::string const &_fname) :
    ice_coupler(_ice_coupler), contract(_contract), fname(_fname)
{
    printf("BEGIN IceWriter::init(%s)\n", fname.c_str());

    ibmisc::Indexing const &indexing(ice_coupler->ice_regridder->gridI->indexing);

    file_initialized = false;

    // Try to be clever about making multi-dimensional arrays
    // in the output according to the grid the user expects.
    dim_names = {"time"};
    counts = {1};
    cur = {0};
    for (size_t i=0; i<indexing.rank(); ++i) {
        // ix goes 0...n-1 for row-major, n-1..0 for col-major
        int ix = indexing.indices()[i];
        dim_names.push_back(string_printf("dim%d", ix));
        counts.push_back(indexing[i].extent);
        cur.push_back(0);
    }

    // Put our output files in this directory, one named per ice sheet.
    boost::filesystem::path ofname(_fname);
    boost::filesystem::create_directory(ofname.parent_path());    // Make sure it exists

//    auto output_dir = boost::filesystem::absolute(
//        boost::filesystem::path(
//            (io == IceCoupler::INPUT ? "ice_model_in" : "ice_model_out")),
//        coupler->gcm_params.run_dir);
//    boost::filesystem::create_directory(output_dir);    // Make sure it exists

    // Set up the output file
    // Create netCDF variables based on details of the coupling contract.xs
    printf("IceWriter opening file %s\n", fname.c_str());

    printf("END IceWriter::init_from_ice_model(%s)\n", fname.c_str());
}

void IceWriter::init_file()
{
    GCMCoupler const *gcm_coupler(ice_coupler->gcm_coupler);
    GCMParams const &gcm_params(gcm_coupler->gcm_params);

    printf("BEGIN IceWriter::init_file(%s)\n", fname.c_str());

    NcIO ncio(fname, NcFile::replace);

    auto one_dims(get_or_add_dims(ncio, {"one"}, {1}));

    // Make time dimension unlimited (-1)
    std::vector<long> dcounts;
    dcounts.reserve(counts.size());
    for (size_t count : counts) dcounts.push_back(count);
    dcounts[0] = -1;        // time needs to be unlimited
    auto dims(get_or_add_dims(ncio, dim_names, dcounts));


//    NcDim one_dim = ncio.nc->addDim("one", 1);
//    NcVar info_var = ncio.nc->addVar("grid", ibmisc::get_nc_type<double>(), one_dims[0]);
    NcVar info_var = get_or_add_var(ncio, "grid", ibmisc::get_nc_type<double>(), one_dims);


    info_var.putAtt("icebin_in", gcm_coupler->icebin_in);
//    info_var.putAtt("variable", coupler->vname);
    info_var.putAtt("ice_sheet", ice_coupler->name());

    NcVar time0_var = get_or_add_var(ncio, "time0", ibmisc::get_nc_type<double>(), one_dims);
    time0_var.putAtt("units", gcm_coupler->time_unit.to_cf());
    time0_var.putAtt("calendar", gcm_coupler->time_unit.cal->to_cf());
    time0_var.putAtt("axis", "T");
    time0_var.putAtt("long_name", "Simulation start time");

    NcVar time_var = get_or_add_var(ncio, "time", ibmisc::get_nc_type<double>(), {dims[0]});
    time_var.putAtt("units", gcm_coupler->time_unit.to_cf());
    time_var.putAtt("calendar", gcm_coupler->time_unit.cal->to_cf());
    time_var.putAtt("axis", "T");
    time_var.putAtt("long_name", "Coupling times");


    for (size_t i=0; i < contract->size(); ++i) {
        VarMeta const &cf = (*contract)[i];
        NcVar var = get_or_add_var(ncio, cf.name, ibmisc::get_nc_type<double>(), dims);
        var.putAtt("units", cf.units);
        var.putAtt("description", cf.description);
    }

    // Put initial time in it...
    time0_var.putVar({0}, {1}, &gcm_coupler->time_start_s);
    ncio.close();

    file_initialized = true;

    printf("END IceWriter::init_file(%s)\n", ice_coupler->name().c_str());
}

template<class T>
T &no_const(const T &val)
{
    return const_cast<T>(val);
}

/** @param index Index of each grid value.
@param vals The values themselves -- could be SMB, Energy, something else...
TODO: More params need to be added.  Time, return values, etc. */
void IceWriter::write(double time_s,
    blitz::Array<double,2> const &valsI)    // valsI[nI, nvars]
{
    GCMCoupler const *gcm_coupler(ice_coupler->gcm_coupler);
    GCMParams const &gcm_params(gcm_coupler->gcm_params);


printf("BEGIN IceWriter::write\n");
    if (!file_initialized) init_file();
    NcIO ncio(fname, NcFile::write);

    // Read index info
    cur[0] = ncio.nc->getDim("time").getSize();

    // Write the current time
    NcVar time_var = ncio.nc->getVar("time");
    time_var.putVar(cur, counts, &time_s);

    // Write the other variables
    size_t nI(valsI.extent(0));
    blitz::Array<double,2> valI_tmp(nI);

    // Sanity check array sizes
    size_t all_count = 1;
    for (auto count : counts) all_count *= count;
    if (all_count != nI) (*icebin_error)(-1,
        "Illegal count to write: %ld vs %ld", all_count, nI);


    for (int ivar=0; ivar < contract->size(); ++ivar) {
        VarMeta const &cf = (*contract)[ivar];

        NcVar ncvar = ncio.nc->getVar(cf.name.c_str());
        for (int i=0; i<valsI.extent(0); ++i) valI_tmp(i) = valsI(i, ivar);
        ncvar.putVar(cur, counts, valI_tmp.data());
    }

    ncio.close();
printf("END IceWriter::write\n");
    
}

}    // namespace icebin

