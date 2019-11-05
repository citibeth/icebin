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
#include <boost/filesystem.hpp>

#include <spsparse/blitz.hpp>
#include <ibmisc/datetime.hpp>
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

namespace icebin {

static double const nan = std::numeric_limits<double>::quiet_NaN();

std::unique_ptr<IceCoupler> new_ice_coupler(NcIO &ncio,
    std::string const &vname, std::string const &sheet_name,
    GCMCoupler const *_gcm_coupler)
{
    std::string vname_sheet(vname + "." + sheet_name);
    auto info_v = get_or_add_var(ncio, vname_sheet + ".info", "int", {});

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
    self->ice_regridder = &*_gcm_coupler->gcm_regridder->ice_regridders().at(sheet_name);
    self->emI_ice.reference(blitz::Array<double,1>(self->ice_regridder->nI()));
    self->emI_ice = nan;
    self->emI_land.reference(blitz::Array<double,1>(self->ice_regridder->nI()));
    self->emI_land = nan;

//    if (rw_full) ncio_blitz(ncio, elevmaskI, true, vname + ".elevmaskI", "double",
//        get_dims(ncio ,{vname + ".gridI.cells.nfull"}));

//    self->ice_constants.init(&_coupler->ut_system);

    self->ncread(ncio, vname_sheet);

    return self;
}

// Default NOP for reconstruct_ice_ivalsI
static void _reconstruct_ice_ivalsI(
    blitz::Array<double,2> &ice_ivalsI,
    double dt) {}

IceCoupler::IceCoupler(IceCoupler::Type _type) :
    type(_type),
    reconstruct_ice_ivalsI(std::bind(
        &_reconstruct_ice_ivalsI, std::placeholders::_1, std::placeholders::_2))
    {}

void IceCoupler::ncread(ibmisc::NcIO &ncio_config, std::string const &vname_sheet)
{
    // General args passed to the ice sheet, regardless of which ice model is being used
    NcVar info_var(ncio_config.nc->getVar(vname_sheet + ".info"));
    get_or_put_att<NcVar,double>(info_var, 'r', "sigma", "double", &sigma[0], 3);
}


IceCoupler::~IceCoupler() {}

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
                (io == 0 ? "icemodel-in.nc" : "icemodel-out.nc"));
            this->writer[io].reset(new IceWriter(
                this, &contract[io], fname.string()));
        }
    }

    // Subclass-specific cold start
    _cold_start(time_base, time_start_s);

    // Allocate
    ice_ovalsI.reference(blitz::Array<double,2>(contract[OUTPUT].size(), nI()));
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

static Eigen::VectorXd to_col_vector(blitz::Array<double,1> const &vec)
{
    Eigen::VectorXd ret(vec.extent(0));
    for (int i=vec.lbound(0), j=0; i <= vec.ubound(0); ++i,++j) ret(j) = vec(i);
    return ret;
}

#if 0
static void mask_result(EigenDenseMatrixT &ret, blitz::Array<double,1> const &wB_b, double fill)
{
    int nB = ret.rows();    // == wB_b.extent(0)
    int nvar = ret.cols();

    // Mask out cells that slipped into the output because they were
    // in the SparseSet; but don't actually get any contribution.
    for (int i=0; i<nB; ++i) {
        if (wB_b(i) != 0.) continue;
printf("wB_b-%d = %f\n", wB_b(i));
        for (int n=0; n<nvar; ++n) ret(i,n) = fill;
    }

}
#endif
// -----------------------------------------------------------
blitz::Array<double,2> IceCoupler::construct_ice_ivalsI(
blitz::Array<double,2> const &gcm_ovalsE0,
std::vector<std::pair<std::string, double>> const &scalars,
double dt,
ibmisc::TmpAlloc &tmp)
{
printf("BEGIN construct_ice_ivalsI(dt=%g)\n", dt);
    auto nE0(gcm_ovalsE0.extent(0));

    // ------------- Form ice_ivalsI
    // Assuming column-major matrices...
    // ice_ivalsI{ik} = IvE0{ij} * (gcm_ovalsE0{jl} * ivg.M{lk} + ivg.b{jk})
    //       Dense    =  Sparse          Dense          Sparse     Repeated
    // *** where:
    // ivg = icei_v_gcmo
    // |i| = # ice grid cells (|I|)
    // |j| = # dense elevation grid cells (|E0|)
    // |k| = # variables in ice_input
    // |l| = # variables in gcm_output

    auto &ice_ivalsI_e(tmp.make<EigenDenseMatrixT>());    // Memory for ice_ivalsI

    // Get the sparse matrix to convert GCM output variables to ice model inputs
    // This will be transposed: M(input, output).  b is a row-vector here.
    auto icei_v_gcmo_T(var_trans_inE.apply_scalars(scalars, 'T'));    // Mxb
//  print_var_trans(icei_v_gcmo_T, var_trans_inE, 'T');

    // Switch from row-major (Blitz++) to col-major (Eigen) indexing
    Eigen::Map<EigenDenseMatrixT> const gcm_ovalsE0_e(
        const_cast<double *>(gcm_ovalsE0.data()),
        gcm_ovalsE0.extent(1), gcm_ovalsE0.extent(0));

    // Ice inputs calculated as the result of a matrix multiplication
    // ice_ivalsI_e is |i| x |k|
    ice_ivalsI_e = (*IvE0->M) * (
        gcm_ovalsE0_e * icei_v_gcmo_T.M + icei_v_gcmo_T.b.replicate(nE0,1) );

    // Alias the Eigen matrix to blitz array
    blitz::Array<double,2> ice_ivalsI(
        ice_ivalsI_e.data(),
        blitz::shape(ice_ivalsI_e.cols(), ice_ivalsI_e.rows()),
        blitz::neverDeleteData);

printf("AA6\n");
    // Continue construction in a contract-specific manner
    reconstruct_ice_ivalsI(ice_ivalsI, dt);

printf("END construct_ice_ivalsI()\n", dt);
    return ice_ivalsI;
}
// -----------------------------------------------------------
/** 
@param do_run True if we are to actually run (otherwise just return ice_ovalsI from current state)
@param gcm_ivalsAE_s Contract inputs for the GCM on the A nad E grid, respectively (1D indexing).
       Accumulates here from many ice shets. _s = sparse indexing
*/
IceCoupler::CoupleOut IceCoupler::couple(
double time_s,
// Values from GCM, passed GCM -> Ice
VectorMultivec const &gcm_ovalsE_s,
// ------- Output Variables (Uses IndexAE::A and IndexAE::E in them)
std::vector<VectorMultivec> &gcm_ivalss_s,                // (accumulate over many ice sheets)
// ------- Flags
bool run_ice)
{
    IceCoupler::CoupleOut ret;

    if (!gcm_coupler->am_i_root()) {
        printf("[noroot] BEGIN IceCoupler::couple(%s)\n", name().c_str());

        // Allocate dummy variables, even though they will only be set on root
        blitz::Array<double,2> ice_ivalsI(contract[INPUT].size(), nI());
        ice_ivalsI = 0;
        ice_ovalsI = 0;
        run_timestep(time_s, ice_ivalsI, ice_ovalsI, run_ice);
        printf("[noroot] END IceCoupler::couple(%s)\n", name().c_str());
        return ret;
    }

    printf("BEGIN IceCoupler::couple(%s)\n", name().c_str());

    // ========== Get Ice Inputs
    // E_s = Elevation grid (sparse indices)
    // E0 = Elevation grid @ beginning of timestep (dense indices)
    // E1 = Elevation grid @ end of timestep (dense indices)
    // All matrices and vectors assumed w/ densified indexing
    // Except _s ending means they use sparse indexing.

    // ------------- Compute dimE0 transformation, if this is the first round
    if (!run_ice) {
        dimE0.reset(new SparseSetT);
        for (size_t i=0; i<gcm_ovalsE_s.size(); ++i) {
            long iE_s(gcm_ovalsE_s.index[i]);
            dimE0->add_dense(iE_s);
        }
    }

    // ------------- Form gcm_ovalsE0
    // Densify gcm_ovalsE_s --> gcm_ovalsE
    // This should ONLY involve iE already mentioned in IvE0;
    // if not, ibmisc_error() will be called inside to_dense()
    blitz::Array<double,2> gcm_ovalsE0(gcm_coupler->gcm_outputsE.size(), dimE0->dense_extent());
    gcm_ovalsE0 = 0;
    for (size_t i=0; i<gcm_ovalsE_s.size(); ++i) {
        long iE_s(gcm_ovalsE_s.index[i]);
        int iE0(dimE0->to_dense(iE_s));   // Can raise error if iE_s not found
        for (int ivar=0; ivar<gcm_ovalsE_s.nvar; ++ivar) {
            gcm_ovalsE0(ivar, iE0) += gcm_ovalsE_s.val(ivar, i);
        }
    }


    // Set up scalars used to instantiate variable conversion matrices
    // NOTE: by_dt=inf on the first call (run_ice=false)
    //       Should be OK because output variables that depend on by_dt
    //       are not needed on the first call.
    double const dt = (time_s - gcm_coupler->last_time_s);
    std::vector<std::pair<std::string, double>> scalars({
        std::make_pair("by_dt", 1.0 / dt)});

    {
        TmpAlloc tmp;    // Allocate variables for the duration of this function
        blitz::Array<double,2> ice_ivalsI(run_ice ?
            construct_ice_ivalsI(gcm_ovalsE0, scalars, dt, tmp) :
            blitz::Array<double,2>(contract[INPUT].size(), nI()));

        // ========= Step the ice model forward
        if (writer[INPUT].get()) {
printf("writing icemodel-in\n");
            writer[INPUT]->write(time_s, ice_ivalsI);
        }
        ice_ovalsI = 0;
        run_timestep(time_s, ice_ivalsI, ice_ovalsI, run_ice);
        if (writer[OUTPUT].get()) {
printf("writing icemodel-out\n");
            writer[OUTPUT]->write(time_s, ice_ovalsI);
        }
    }

    // ========== Update regridding matrices
    int emI_ice_ix = standard_names[OUTPUT].at("elevmask_ice");
    blitz::Array<double,1> out_emI_ice(ice_ovalsI(emI_ice_ix, blitz::Range::all()));

    int emI_land_ix = standard_names[OUTPUT].at("elevmask_land");
    blitz::Array<double,1> out_emI_land(ice_ovalsI(emI_land_ix, blitz::Range::all()));

    // Check that elevmaskI is an alias for variable #elevmaskI_ix in ice_ovalsI
    if (&ice_ovalsI(emI_ice_ix,0) != &out_emI_ice(0)) (*icebin_error)(-1,
        "ice_ovalsI <%p> != emI_ice <%p>\n", &ice_ovalsI(emI_ice_ix,0), &out_emI_ice(0));
    if (&ice_ovalsI(emI_land_ix,0) != &out_emI_land(0)) (*icebin_error)(-1,
        "ice_ovalsI <%p> != emI_land <%p>\n", &ice_ovalsI(emI_land_ix,0), &out_emI_land(0));

    emI_ice = out_emI_ice;    // Copy
    emI_land = out_emI_land;    // Copy
    GCMRegridder *gcmr(&*gcm_coupler->gcm_regridder);
    int sheet_index = gcmr->ice_regridders().index.at(name());
    std::unique_ptr<RegridMatrices_Dynamic> rm(gcmr->regrid_matrices(sheet_index, emI_ice));

    // ------ Update E1vE0 translation between old and new elevation classes
    //        (global for all ice sheets)
    // A SparseSet that is identity for the entire range of I
    SparseSetT dimI(id_sparse_set<SparseSetT>(nI()));

    // _nc means "No Correct" for changes in area due to projections
    // See commit d038e5cb for deeper explanation
    std::unique_ptr<SparseSetT> dimE1(new SparseSetT);
    ret.E1vI_unscaled_nc = rm->matrix_d("EvI", {&*dimE1, &dimI},
        RegridParams(false, false, {0,0,0}));    // scale=f, correctA=f
printf("dimE1 extents 1: %p %ld %ld\n", &*dimE1, (long)dimE1->sparse_extent(), (long)dimE1->dense_extent());

    // ========= Compute gcm_ivalsE
    SparseSetT dimA1;
    auto A1vI_unscaled(rm->matrix_d("AvI", {&dimA1, &dimI},
        RegridParams(false, true, {0,0,0})));    // scale=f, correctA=t

    // Do it once for _E variables and once for _A variables.
    std::vector<linear::Weighted_Eigen *> AE1vIs(gcm_ivalss_s.size());
        AE1vIs[(int)IndexAE::A] = &*A1vI_unscaled;
        AE1vIs[(int)IndexAE::E] = &*ret.E1vI_unscaled_nc;

    for (int iAE=(int)IndexAE::A; iAE <= (int)IndexAE::E; ++iAE) {

        // Assuming column-major matrices...
        // gcm_ivalsX{jn} = X1vI{ji} * (ice_ovalsI{im} * gvi.M{mn} + gvi.b{in})
        //      Dense        Sparse         Dense          Sparse     Repeated

        // *** where:
        // gvi = gcmi_v_iceo (unit/variable conversion)
        // |i| = # ice grid cells (|I|)
        // |j| = # used elevation/atmosphere grid cells (|X|)
        // |m| = # variables in ice_output
        // |n| = # variables in gcm_input

        // Get the sparse matrix to convert ice model output variables to GCM inputs
        // This will be transposed: M(input, output).  b is a row-vector here.
        auto gcmi_v_iceo_T(var_trans_outAE[iAE].apply_scalars(scalars, 'T'));
//        print_var_trans(gcmi_v_iceo_T, var_trans_outAE[iAE], 'T');

        // Switch from row-major (Blitz++) to col-major (Eigen) indexing
        Eigen::Map<EigenDenseMatrixT> ice_ovalsI_e(
            ice_ovalsI.data(), ice_ovalsI.extent(1), ice_ovalsI.extent(0));

        // ----------- Sanity check: There should not be any NaNs...
        bool hasnan = false;
        for (auto ii(begin(*AE1vIs[iAE]->M)); ii != end(*AE1vIs[iAE]->M); ++ii) {
            if (std::isnan(ii->value())) {
            	    printf("nan found: AvI[%d, %d] = %g\n", ii->row(), ii->col(), ii->value());
                hasnan = true;
            }
        }
        for (int j=0; j<ice_ovalsI_e.cols(); ++j) {
            if (contract[OUTPUT][j].flags & contracts::ALLOW_NAN) continue;
            for (int i=0; i<ice_ovalsI_e.rows(); ++i) {
                auto val(ice_ovalsI_e(i,j));
                if (std::isnan(val)) {
                    printf("nan found: ice_ovalsI_e[%d, %d] = %g\n", i, j, val);
                    hasnan = true;
                }
            }
        }
        for (auto ii(begin(gcmi_v_iceo_T.M)); ii != end(gcmi_v_iceo_T.M); ++ii) {
            if (std::isnan(ii->value())) {
                printf("nan found: gcmi_v_iceo_T.M[%d, %d] = %g\n", ii->row(), ii->col(), ii->value());
                hasnan = true;
            }
        }

        for (int j=0; j<gcmi_v_iceo_T.b.cols(); ++j) {
        for (int i=0; i<gcmi_v_iceo_T.b.rows(); ++i) {
            auto val(gcmi_v_iceo_T.b(i,j));
            if (std::isnan(val)) {
                printf("nan found: gcmi_v_iceo_T.b[%d, %d] = %g\n", i, j, val);
                hasnan = true;
            }
        }}

        if (hasnan) (*icebin_error)(-1, "At least one NaN detected!");
        // -------------------------- END Sanity Check

        // Regrid while recombining variables
        // (Do not need to use Weighted_Eigen::apply(), since this is not IvE)
        EigenDenseMatrixT gcm_ivalsX((*AE1vIs[iAE]->M) * (
            ice_ovalsI_e * gcmi_v_iceo_T.M + gcmi_v_iceo_T.b.replicate(nI(),1) ));
        // Sparsify while appending to the global VectorMultivec
        // (Transposes order in memory)
        std::vector<double> vals(gcm_ivalss_s[iAE].nvar);
        for (int jj=0; jj < gcm_ivalsX.rows(); ++jj) {
            auto jj_s(AE1vIs[iAE]->dims[0]->to_sparse(jj));
            int nn = 0;
            for (; nn < gcm_ivalsX.cols(); ++nn) {
                vals[nn] = gcm_ivalsX(jj,nn);
            }
            gcm_ivalss_s[iAE].add(jj_s, vals, AE1vIs[iAE]->wM(jj));

        }
    }        // iAE
    // Compute IvE (for next timestep)
    std::unique_ptr<linear::Weighted_Eigen> IvE1(
        rm->matrix_d("IvE", {&dimI, &*dimE1},
        RegridParams(true, true, {0,0,0}))); // scale=t, correctA=t

    // wIvE0.reference(IvE1->wM);


    // Record our matrices for posterity
    {auto fname(
        boost::filesystem::path(output_dir) / 
        ("regrids-" + ice_regridder->name() + "-" + gcm_coupler->sdate(time_s) + ".nc"));
    NcIO ncio(fname.string(), NcFile::replace);

        // Write matrices as their dense subspace versions, not the sparsified versions.
        dimI.ncio(ncio, "dimI");
        dimA1.ncio(ncio, "dimA");
        dimE1->ncio(ncio, "dimE");

        AE1vIs[GridAE::E]->ncio(ncio, "EvI_unscaled_nc", {"dimE", "dimI"});
        AE1vIs[GridAE::A]->ncio(ncio, "AvI_unscaled", {"dimA", "dimI"});
        IvE1->ncio(ncio, "IvE", {"dimI", "dimE"});
    }

    // Save stuff for next time around
    if (run_ice) {    // But not if first time
        ret.dimE0 = std::move(this->dimE0);
        ret.IvE0 = std::move(IvE0->M);
    }
    this->dimE0 = std::move(dimE1);
    this->IvE0 = std::move(IvE1);

    printf("END IceCoupler::couple(%s)\n", name().c_str());
    return ret;
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

    ibmisc::Indexing const &indexing(ice_coupler->ice_regridder->agridI.indexing);

    file_initialized = false;

    // Try to be clever about making multi-dimensional arrays
    // in the output according to the grid the user expects.
    dim_names = {"time"};
    counts = {1};
    cur = {0};
    for (size_t i=0; i<indexing.rank(); ++i) {

        // ix goes 0...n-1 for row-major, n-1..0 for col-major
        int ix = indexing.indices()[i];
        dim_names.push_back(indexing[ix].name);
        counts.push_back(indexing[ix].extent);
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
    auto &time_unit(gcm_coupler->time_unit);

    printf("BEGIN IceWriter::init_file(%s)\n", fname.c_str());

    NcIO ncio(fname, NcFile::replace);

    auto one_dims(get_or_add_dims(ncio, {"one"}, {1}));

    // Make time dimension unlimited (-1)
    std::vector<long> dcounts;
    dcounts.reserve(counts.size());
    for (size_t count : counts) dcounts.push_back(count);
    dcounts[0] = -1;        // time needs to be unlimited
    auto dims(get_or_add_dims(ncio, dim_names, dcounts));


    NcVar info_var = get_or_add_var(ncio, "grid", ibmisc::get_nc_type<double>(), one_dims);

    info_var.putAtt("icebin_in", gcm_coupler->icebin_in);
    info_var.putAtt("ice_sheet", ice_coupler->name());

    NcVar time0_var = get_or_add_var(ncio, "time0", ibmisc::get_nc_type<double>(), one_dims);
    time0_var.putAtt("units", time_unit.to_cf());
    time0_var.putAtt("calendar", time_unit.cal->to_cf());
    time0_var.putAtt("axis", "T");
    time0_var.putAtt("long_name", "Simulation start time");

    NcVar time0_txt = get_or_add_var(ncio, "time0.txt", "char",
        get_or_add_dims(ncio, {"iso8601.length"}, {iso8601_length}));

    NcVar time_var = get_or_add_var(ncio, "time", ibmisc::get_nc_type<double>(), {dims[0]});
    time_var.putAtt("units", time_unit.to_cf());
    time_var.putAtt("calendar", time_unit.cal->to_cf());
    time_var.putAtt("axis", "T");
    time_var.putAtt("long_name", "Coupling times");

    NcVar time_txt = get_or_add_var(ncio, "time.txt", "char",
        get_or_add_dims(ncio, {dim_names[0], "iso8601.length"}, {-1, iso8601_length}));


    for (size_t i=0; i < contract->size(); ++i) {
        VarMeta const &cf = (*contract)[i];
        NcVar var = get_or_add_var(ncio, cf.name, ibmisc::get_nc_type<double>(), dims);
        var.putAtt("units", cf.units);
        var.putAtt("description", cf.description);
    }

    // Put initial time in it...
    double const &time_start_s(gcm_coupler->time_start_s);
    time0_var.putVar({0}, {1}, &time_start_s);
    time0_txt.putVar({0}, {iso8601_length},
        to_iso8601(time_unit.to_datetime(time_start_s)).c_str());

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
    blitz::Array<double,2> const &valsI)    // valsI[nvars, nI]
{
    GCMCoupler const *gcm_coupler(ice_coupler->gcm_coupler);
    GCMParams const &gcm_params(gcm_coupler->gcm_params);
    auto &time_unit(gcm_coupler->time_unit);

    if (!file_initialized) init_file();
    NcIO ncio(fname, NcFile::write);

    // Read index info
    cur[0] = ncio.nc->getDim("time").getSize();

    // Write the current time
    NcVar time_var = ncio.nc->getVar("time");
    time_var.putVar(cur, counts, &time_s);

    NcVar time_txt = ncio.nc->getVar("time.txt");
    time_txt.putVar({cur[0],0}, {counts[0],iso8601_length},
        to_iso8601(time_unit.to_datetime(time_s)).c_str());

    // Write the other variables
    size_t nI(valsI.extent(1));
    blitz::Array<double,1> valI_tmp(nI);

    // Sanity check array sizes
    size_t all_count = 1;
    for (auto count : counts) all_count *= count;
    if (all_count != nI) (*icebin_error)(-1,
        "Illegal count to write: %ld vs %ld", all_count, nI);

    for (int ivar=0; ivar < contract->size(); ++ivar) {
        VarMeta const &cf = (*contract)[ivar];

        NcVar ncvar = ncio.nc->getVar(cf.name.c_str());
        for (int i=0; i<valsI.extent(1); ++i) valI_tmp(i) = valsI(ivar, i);
        ncvar.putVar(cur, counts, valI_tmp.data());
    }

    ncio.close();
}

}    // namespace icebin

