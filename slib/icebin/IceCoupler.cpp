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
#include <spsparse/multiply_sparse.hpp>
#include <spsparse/sort.hpp>
#include <spsparse/eigen.hpp>

#ifdef USE_PISM
#include <icebin/pism/IceCoupler_PISM.hpp>
#endif

using namespace spsparse;
using namespace ibmisc;
using namespace netCDF;
//using namespace std;


namespace icebin {

std::unique_ptr<IceCoupler> new_ice_coupler(NcIO &ncio, std::string vname,
    GCMCoupler const *_gcm_coupler, IceRegridder *_ice_regridder)
{
    std::string vn(vname + ".info");
    auto info_v = get_or_add_var(ncio, vn, "int64", {});

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
    self->gcm_coupler = _gcm_coupler;
    self->ice_regridder = _ice_regridder;
//    self->ice_constants.init(&_coupler->ut_system);

    self->ncread(ncio, vname);

    return self;
}


IceCoupler::~IceCoupler() {}
// ==========================================================
static double const nan = std::numeric_limits<double>::quiet_NaN();

// ------------------------------------------------------------

/** Returns just the matrix part of regridding.  Also sets scale and
    correctA.  Local convenience function. */
inline SparseMatrix regrid_M(
RegridMatrices &rm,
std::string const &spec_name)
{
    std::unique_ptr<WeightedSparse> WS(rm.regrid(spec_name, true, true));
    EigenSparseMatrix &Meig(*WS->M);
    SparseMatrix Mout;
    spcopy(Mout, Meig, true);
    return std::move(Mout);
}



// ==============================================================
void IceCoupler::cold_start(
        ibmisc::Datetime const &time_base,
        double time_start_s)
{
    // Set up writers
    for (int io=0; io<2; ++io) {
        this->writer[io].reset(new IceWriter(
            this, &contract[io],
            name() + (io == 0 ? "_in.nc" : "_out.nc")));
    }

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
VectorMultivec const &gcm_ovalsE,
GCMInput &out,    // Accumulate matrices here...
bool run_ice,
bool am_i_root)
{
    if (!am_i_root) {
        // Allocate dummy variables, even though they will only be set on root
        blitz::Array<double,2> ice_ivalsI(ndata(), contract[INPUT].size());
        blitz::Array<double,2> ice_ovalsI(ndata(), contract[OUTPUT].size());
        ice_ivalsI = 0;
        ice_ovalsI = 0;
        run_timestep(time_s, ice_ivalsI, ice_ovalsI, run_ice, am_i_root);
        return;
    }

    // Store regridding matrices for the last timestep, which we will
    // need to create ice_ivals
//    std::unique_ptr<WeightedSparse> IvEd0(std::move(IvEd));    // Dense A and E
//    std::unique_ptr<WeightedSparse> IvAd0(std::move(IvAd));

    // ========== Get Ice Inputs
    blitz::Array<double,2> ice_ivalsI(ndata(), contract[INPUT].size());
    blitz::Array<double,2> ice_ovalsI(ndata(), contract[OUTPUT].size());
    ice_ivalsI = 0;
    ice_ovalsI = 0;

    // Save old IvE for now
    auto IvE0(std::move(IvE));

    // ------------- Form ice_ivalsI
    // Densify gcm_ovals and IvE
    SparseMatrix IvE0_M; spcopy(IvE0_M, *IvE0->M);
    SparseSet dimE;
    SparseMatrix IvEd0;
    dimE.add_sorted(IvE0_M.dim_begin(1), IvE0_M.dim_end(1));
//    densify_one_dim(IvEd0, IvE0, dimE, 1);
//    sparse_copy(IvEd0, IvE0,
//        SparseTransform::TO_DENSE,
//        make_array((decltype(dimE) *)0, &dimE));

    // Densify gcm_ovalsE --> gcm_ovalsEd
    // This should ONLY involve iE already mentioned in IvE0;
    // if not, there will be an exception.
    blitz::Array<double,2> gcm_ovalsEd(dimE.sparse_extent(), gcm_coupler->gcm_outputsE.size());
    gcm_ovalsEd = 0;
    for (size_t i=0; i<gcm_ovalsE.index.size(); ++i) {
        auto iE(gcm_ovalsE.index[i]);
        int iEd (dimE.to_dense(iE));
        for (int ivar=0; ivar<gcm_ovalsE.size(); ++i) {
            gcm_ovalsEd(iEd, ivar) += gcm_ovalsE.val(ivar, i);
        }
    }

    // Get the CSR sparse matrix to convert GCM outputs to ice model inputs
    

    std::vector<std::pair<std::string, double>> scalars({
        std::make_pair("by_dt", 1.0 / (time_s - gcm_coupler->last_time_s)),
        std::make_pair("unit", 1.0)});
    CSRAndUnits icei_v_gcmo(var_trans_inE.apply_scalars(scalars));

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
    for (auto iIvEd0 = IvEd0.begin(); iIvEd0 != IvEd0.end(); ++iIvEd0) {
        int ii(iIvEd0.index(0));
        int jj(iIvEd0.index(1));
        auto IvEd0_ij(iIvEd0.val());

        // Transform units on the input while multiplying by M
        for (int kk = 0; kk < contract[INPUT].size(); ++kk) {
            double ice_ivalsEd_jk = 0;
            std::vector<std::pair<int, double>> const &row(icei_v_gcmo.mat[kk]);
            for (auto rowk_iter=row.begin(); rowk_iter != row.end(); ++rowk_iter) {
                int ll(rowk_iter->first);
                double icei_v_gcmo_kl(rowk_iter->second);
                ice_ivalsEd_jk += icei_v_gcmo_kl * gcm_ovalsEd(jj, ll);
            }
            ice_ivalsI(ii, kk) += IvEd0_ij * ice_ivalsEd_jk;
        }
    }

    // ========= Step the ice model forward
    if (writer[INPUT].get()) writer[INPUT]->write(time_s, ice_ivalsI);
    run_timestep(time_s, ice_ivalsI, ice_ovalsI, run_ice, am_i_root);
    if (writer[OUTPUT].get()) writer[OUTPUT]->write(time_s, ice_ovalsI);

    // ========== Update regridding matrices
    auto elevI(get_elevI());    // blitz
    ice_regridder->set_elevI(elevI);
    RegridMatrices rm(ice_regridder);

    // Compute IvE (for next timestep)
    auto IvE(regrid_M(rm, "IvE"));
    dimE.add_sorted(IvE.dim_begin(1), IvE.dim_end(1));

    // Compute regrid matrices we need now
    auto EvI(rm.regrid("EvI", true, true));
    auto AvI(rm.regrid("AvI", true, true));
    auto AvE(rm.regrid("AvE", true, true));
    SparseMatrix &out_AvE1(out.AvE1);
    sparse_copy(out_AvE1, *AvE->M,
        SparseTransform::TO_SPARSE,
        make_array(&AvE->dims[0], &AvE->dims[1]));

    SparseVector &out_wAvE1(out.wAvE1);    // Get around bug in gcc@4.9.3
    sparse_copy(out_wAvE1, AvE->weight,
        SparseTransform::TO_SPARSE,
        make_array(&AvE->dims[0]), true);


    // Accumulate global things for all ice sheets....
    SparseMatrix &out_E1vE0(out.E1vE0);
    Eigen::SparseMatrix<double> _prod1(*EvI->M * *IvE0->M);
    sparse_copy(out_E1vE0, _prod1,
        SparseTransform::TO_SPARSE,
        make_array(&EvI->dims[0], &dimE0));   // copy eigen matrix

#if 0
    Eigen::SparseMatrix<double> M(*EvI->M);
    Eigen::VectorXd xx(to_col_vector(elevI));
    Eigen::VectorXd yy(M * xx);

    Eigen::VectorXi xxx(*EvI->M * to_col_vector(elevI));
#endif

    SparseVector &out_elevE1(out.elevE1);
    Eigen::VectorXd _prod2(*EvI->M * to_col_vector(elevI));
    sparse_copy(out_elevE1, _prod2,
        SparseTransform::TO_SPARSE,
        make_array(&EvI->dims[0]));   // copy eigen column vector

    // ========= Compute gcm_ivalsE = EvI * vt * ice_ovals

    std::array<WeightedSparse *, GridAE::count> XvIs;
    XvIs[GridAE::E] = &*EvI;
    XvIs[GridAE::A] = &*AvI;
    std::vector<double> vals(contract.size());

    // Do it once for _E variables and once for _A variables.
    for (int iAE=0; iAE < GridAE::count; ++iAE) {
        CSRAndUnits gcmi_v_iceo(var_trans_outAE[iAE].apply_scalars(scalars));

        VarSet const &contract(gcm_coupler->gcm_inputsAE[iAE]);
        WeightedSparse &XvI_ws(*XvIs[iAE]);
        SparseMatrix XvI;
        sparse_copy(XvI, *XvI_ws.M,
            SparseTransform::TO_SPARSE,
            make_array(&XvI_ws.dims[0], &XvI_ws.dims[1]));
        VectorMultivec &gcm_ivalsX(out.gcm_ivalsAE[iAE]);

        // gcm_ivalsX_{jn} = XvI_{ji} * gcmi_v_iceo_{nm} * ice_ovalsI_{im}
        //       where:
        // |i| = # ice grid cells (|I|)
        // |j| = # elevation/atmosphere grid cells (|X|)
        // |m| = # variables in ice_output
        // |n| = # variables in gcm_input

        // Do the multiplication
        for (auto iXvI(XvI.begin()); iXvI != XvI.end(); ++iXvI) {
            auto jj(iXvI.index(0));
            auto ii(iXvI.index(1));
            auto XvI_ji(iXvI.val());

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
                vals[nn] = XvI_ji * gcm_ivalsX_in;
            }
            gcm_ivalsX.add(jj, vals);
        }
    }

    dimE0 = std::move(dimE);
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
    strides = {1};
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
    auto dims(get_or_add_dims(ncio, dim_names, counts));


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
    strides[0] = &valsI(1,0) - &valsI(0,0);
    for (int ivar=0; ivar < contract->size(); ++ivar) {
        VarMeta const &cf = (*contract)[ivar];

        NcVar ncvar = ncio.nc->getVar(cf.name.c_str());
        double const *data(&valsI(0,ivar));
        ncvar.putVar(cur, counts, strides, data);
    }

    ncio.close();
printf("END IceWriter::write\n");
    
}

}    // namespace icebin

