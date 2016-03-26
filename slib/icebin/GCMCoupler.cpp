/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <mpi.h>        // Intel MPI wants to be first
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <icebin/GCMCoupler.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/contracts/contracts.hpp>
#include <spsparse/multiply_sparse.hpp>

#ifdef USE_PISM
#include <icebin/pism/IceModel_PISM.hpp>
#endif

using namespace spsparse;
using namespace ibmisc;
using namespace netCDF;

namespace icebin {

std::unique_ptr<IceModel> new_ice_model(IceModel::Type type,
    GCMCoupler const *_coupler, IceRegridder *_sheet)
{
    std::unique_ptr<IceModel> ice_model;

    switch(type.index()) {
#if 0
        case IceModel::Type::DISMAL :
            ice_model.reset(new IceModel_DISMAL);
        break;
#endif
        case IceModel::Type::WRITER :
            ice_model.reset(new IceModel_Writer);
        break;
#ifdef USE_PISM
        case IceModel::Type::PISM :
            ice_model.reset(new gpism::IceModel_PISM);
        break;
#endif
        default :
            (*icebin_error)(-1,
                "Unknown IceModel::Type %s", type.str());
    }


    // Do basic initialization...
    ice_model->coupler = _coupler;
    ice_model->sheet = _sheet;
    ice_model->ice_constants.init(&_coupler->ut_system);

    return ice_model;


}

std::unique_ptr<IceModel> new_ice_model(NcIO &ncio, std::string vname,
    GCMCoupler const *_coupler, IceRegridder *_sheet)
{
    std::string vn(vname + ".info");
    auto info_v = get_or_add_var(ncio, vn, netCDF::ncInt64, {});

    IceModel::Type type;
    get_or_put_att_enum(info_v, ncio.rw, "ice_model", type);

    return new_ice_model(type, _coupler, _sheet);
}

void IceModel::ncread(
ibmisc::NcIO &ncio, std::string const &vname_sheet)
{
    gcm_per_ice_sheet_params = coupler->read_gcm_per_ice_sheet_params(ncio, vname_sheet);
}



IceModel::~IceModel() {}

#if 0
VarSet *IceModel::new_VarSet() {
    _extra_contracts.push_back(
        std::unique_ptr<VarSet>(
        new VarSet()));
    return _extra_contracts[_extra_contracts.size()-1].get();
}
#endif

// ==========================================================
bool IceModel::am_i_root() const
    { return coupler->am_i_root(); }

/** Allocate vectors in preparation of calling an ice model. */
void IceModel::allocate_ice_ovals_I()
{
    // Check for program errors
    if (!coupler->am_i_root()) (*icebin_error)(-1,
        "IceModel::allocate_ice_ovals_I() should only be called from GCM root MPI node.  Fix the code.");

    if (ice_ovals_I.size() != 0) (*icebin_error)(-1,
        "[%d] IceModel::allocate_ice_ovals_I(): called twice without a free() inbetween.  Fix the code. (old size is %ld)\n", coupler->gcm_params.gcm_rank, ice_ovals_I.size());

    // Allocate for direct output from ice model
    VarSet const &ocontract(contract[IceModel::OUTPUT]);
    for (int i=0; i < ocontract.size(); ++i)
        ice_ovals_I.push_back(blitz::Array<double,1>(ndata()));
}


/** Allocate in preparation of var transformations (but not regridding yet) */
void IceModel::allocate_gcm_ivals_I()
{
    // Check for program errors
    if (!coupler->am_i_root()) (*icebin_error)(-1,
        "IceModel::allocate_ice_ivals_I() should only be called from GCM root MPI node.  Fix the code.\n");

    if (gcm_ivals_I.size() != 0) (*icebin_error)(-1,
        "IceModel::allocate_gcm_ivals_I(): called twice without a free() inbetween.  Fix the code.\n");


    VarSet const &gcm_inputs(coupler->gcm_inputs);
    for (int i=0; i < gcm_inputs.size(); ++i) {
        gcm_ivals_I.push_back(blitz::Array<double,1>(ndata()));
    }
}

/** Free portions not needed after finished calling ice model and
applying variable transform.  This will be variables desired on
anything other than the ELEVATION grid. */
void IceModel::free_ice_ovals_I()
{
    // Check for program errors
    if (!coupler->am_i_root()) (*icebin_error)(-1,
        "IceModel::free_ice_ovals_I() should only be called from GCM root MPI node.  Fix the code.\n");

    ice_ovals_I.clear();
}

/** Free all memory used by this.  Called when we're done with a coupling timestep. */
void IceModel::free_ovals_ivals_I()
{
    // Check for program errors
    if (!coupler->am_i_root()) (*icebin_error)(-1,
        "IceModel::free_ovals_ovals_I() should only be called from GCM root MPI node.  Fix the code.\n");

    ice_ovals_I.clear();
    gcm_ivals_I.clear();
}

// -----------------------------------------------------------
/** Allocates and sets gcm_ivals_I variable
@param mask Control which GCM input variables to set (according to, etc, INITIAL flag in contract) */
void IceModel::set_gcm_inputs(unsigned int mask)
{
  printf("BEGIN IceModel::set_gcm_inputs()\n");
    allocate_gcm_ivals_I();

    // Compute the variable transformation
    ibmisc::VarTransformer &vt(var_transformer[IceModel::OUTPUT]);
    ibmisc::CSRAndUnits trans = vt.apply_scalars({
        std::make_pair("unit", 1.0)});

    VarSet const &gcm_inputs(coupler->gcm_inputs);

    // Apply the variable transformation
    for (int xi=0; xi<vt.dim(ibmisc::VarTransformer::OUTPUTS).size(); ++xi) {   // xi is index of output variable
        gcm_ivals_I[xi] = 0;    // Vector operation: clear before sum
        VarMeta const &cf(gcm_inputs[xi]);

        if ((cf.flags & mask) != mask) continue;

        // Consider each output variable separately...
        std::vector<std::pair<int, double>> const &row(trans.mat[xi]);
        for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
            int xj = xjj->first;        // Index of input variable
            double io_val = xjj->second;    // Amount to multiply it by
            gcm_ivals_I[xi] += ice_ovals_I[xj] * io_val;        // blitz++ vector operation
        }
        gcm_ivals_I[xi] += trans.units[xi];
    }

    printf("END IceModel::set_gcm_inputs()\n");
}

// -----------------------------------------------------------









// ==========================================================
/** @param nc The Glint22 configuration file */
void GCMCoupler::ncread(
    std::string const &_fname,
    std::string const &_vname,
    ibmisc::Domain<int> &&_domainA)
{
    printf("BEGIN GCMCoupler::read_from_netcdf() %s\n", vname.c_str());

    fname = _fname;
    vname = _vname;
    domainA = std::move(_domainA);

    NcIO ncio(fname, NcFile::read);

    // Load the MatrixMaker (filtering by our domain, of course)
    // Also load the ice sheets
    regridder.ncio(ncio, vname);
    regridder.filter_cellsA(domainA);

    // The root node needs to make the full regridding matrices, not just those
    // for its own domain.  Therefore, create a second MatrixMaker for it.
    if (am_i_root()) {
        regridder_full.reset(new GCMRegridder);
        regridder_full->ncio(ncio, vname);
    }

    // Read gcm_out_file, an optional variable telling the GCM-specific
    // part of IceBin to write out exactly what it sees coming from the GCM
    // (so it can be replayed later with desm)
    auto info_v = get_or_add_var(ncio, vname + ".info", netCDF::ncInt64, {});

    gcm_out_file = "";
    get_or_put_att(info_v, ncio.rw, "gcm_out_file", gcm_out_file, false);
    if (gcm_out_file.length() > 0) {
        gcm_out_file = boost::filesystem::absolute(
            boost::filesystem::path(gcm_out_file),
            gcm_params.run_dir).string();
    }

    gcm_in_file = "";
    get_or_put_att(info_v, ncio.rw, "gcm_in_file", gcm_in_file, false);
    if (gcm_in_file.length() > 0) {
        gcm_in_file = boost::filesystem::absolute(
            boost::filesystem::path(gcm_in_file),
            gcm_params.run_dir).string();
    }



#if 1
    std::cout << "========= GCM Constants" << std::endl;
    std::cout << gcm_constants;
    std::cout << "========= GCM Outputs" << std::endl;
    std::cout << gcm_outputs;
    std::cout << "========= GCM Inputs" << std::endl;
    std::cout << gcm_inputs;


#endif

    ice_models.clear();
    for (size_t i=0; i < regridder.sheets.size(); ++i) {
        IceRegridder *sheet = &*regridder.sheets[i];
        std::string vname_sheet = vname + "." + sheet->name();

        // Create an IceModel corresponding to this IceSheet.
        std::unique_ptr<IceModel> ice_model(new_ice_model(ncio, vname_sheet, this, sheet));
        ice_model->ncread(ncio, vname_sheet);

        contracts::setup(*this, *ice_model);    // where does this go w.r.t ncread() and upate?
        ice_model->update_ice_sheet(ncio, vname_sheet);

        // ---- Create the affiliated writers
        // Writer for ice model input
        auto iwriter(dynamic_cast_unique_ptr<IceModel_Writer, IceModel>(
            new_ice_model(IceModel::Type::WRITER, this, sheet)));
        iwriter->init(IceModel::INPUT, &*ice_model);        // Writer-specific initialization
        ice_model->_iwriter = std::move(iwriter);

        // Writer for ice model output
        auto owriter(dynamic_cast_unique_ptr<IceModel_Writer, IceModel>(
            new_ice_model(IceModel::Type::WRITER, this, sheet)));
        owriter->init(IceModel::OUTPUT, &*ice_model);       // Writer-specific initialization
        ice_model->_owriter = std::move(owriter);

        ice_models.push_back(std::move(ice_model));
    }

    ncio.close();
    printf("END GCMCoupler::read_from_netcdf()\n");
}



void GCMCoupler::set_start_time(
    ibmisc::time::tm const &time_base,
    double time_start_s)
{
    gcm_params.set_start_time(time_base, time_start_s);

    for (size_t sheetix=0; sheetix < ice_models.size(); ++sheetix) {
        auto &ice_model(ice_models[sheetix]);

        ice_model->start_time_set();

        // Handle writing stuff, if the
        IceModel_Writer *iwriter = ice_model->iwriter();    // The affiliated input-writer (if it exists).
        if (iwriter) iwriter->start_time_set();

        IceModel_Writer *owriter = ice_model->owriter();
        if (owriter) owriter->start_time_set();

        // This function is the last phase of initialization.
        // Only now can we assume that the contracts are fully set up.
#if 1
        // Print out the contract and var transformations
        std::cout << "========= Contract for " << ice_model->name() << std::endl;
        std::cout << "---- GCM->Ice     Output Variables:" << std::endl;
        std::cout << ice_model->contract[IceModel::INPUT];
        std::cout << "TRANSFORMATIONS:" << std::endl;
        std::cout << ice_model->var_transformer[IceModel::INPUT];
        std::cout << "---- Ice->GCM     Output Variables:" << std::endl;
        std::cout << ice_model->contract[IceModel::OUTPUT];
        std::cout << "TRANSFORMATIONS:" << std::endl;
        std::cout << ice_model->var_transformer[IceModel::OUTPUT];
#endif
    }

}



// ===================================================
// SMBMsg

MPI_Datatype SMBMsg::new_MPI_struct(int nfields)
{
    int nele = 2 + nfields;
    int blocklengths[] = {1, 1, nfields};
    MPI_Aint displacements[] = {offsetof(SMBMsg,sheetix), offsetof(SMBMsg,iI), offsetof(SMBMsg, vals)};
    MPI_Datatype types[] = {MPI_INT, MPI_INT, MPI_DOUBLE};
    MPI_Datatype ret;
    MPI_Type_create_struct(3, blocklengths, displacements, types, &ret);
    MPI_Type_commit(&ret);
    return ret;
}

/** for use with qsort */
int SMBMsg::compar(void const *a, void const *b)
{
    SMBMsg const *aa = reinterpret_cast<SMBMsg const *>(a);
    SMBMsg const *bb = reinterpret_cast<SMBMsg const *>(b);
    int cmp = aa->sheetix - bb->sheetix;
    if (cmp != 0) return cmp;
    return aa->iI - bb->iI;
}

// ===================================================
// GCMCoupler

/** Parameters to the call_ice_model() method.  Calls to the
ice model are planned first, then executed separately. */
class CallIceModelParams {
public:
    int sheetix;
    SMBMsg *begin;
    SMBMsg *next;

    CallIceModelParams() {}

    CallIceModelParams(int _sheetix, SMBMsg *_begin, SMBMsg *_next) :
        sheetix(_sheetix), begin(_begin), next(_next) {}
};

/** PROTECTED method.  Called by all MPI nodes. */
void GCMCoupler::call_ice_model(
    IceModel *model,
    int sheetix,
    double time_s,
    ibmisc::DynArray<SMBMsg> &rbuf,
    SMBMsg *begin, SMBMsg *end)     // Messages have the MAX number of fields for any ice model contract
{
    // ----------------- Construct input arrays
    // The number of fields for THIS ice sheet will be <= the number
    // of fields in SMBMsg
    int nfields = model->contract[IceModel::INPUT].size();

printf("BEGIN GCMCoupler::call_ice_model(%s, nfields=%ld)\n", model->name().c_str(), nfields);

    // Construct indices vector
    blitz::TinyVector<int,1> shape(rbuf.diff(end, begin));
    blitz::TinyVector<int,1> stride(rbuf.ele_size / sizeof(int));
    blitz::Array<int,1> indices(&begin->iI,
        shape, stride, blitz::neverDeleteData);

    // Construct input values vectors: sparse vectors pointing into the SMBMsgs
    std::vector<blitz::Array<double,1>> ivals2;
    stride[0] = rbuf.ele_size / sizeof(double);
    for (int i=0; i<nfields; ++i) {
        SMBMsg &rbegin(*begin);     // Reference to begin

        ivals2.push_back(blitz::Array<double,1>(
            &rbegin[i], shape, stride, blitz::neverDeleteData));
    }


    // -------------- Run the model
    // Record exactly the same inputs that this ice model is seeing.
    IceModel_Writer *iwriter = model->iwriter();
    if (iwriter) iwriter->run_timestep(time_s, indices, ivals2);

    // Now call to the ice model
//  printf("[%d] BEGIN GCMCoupler::call_ice_model() calling model->run_timestep()\n", rank);
    model->run_timestep(time_s, indices, ivals2);
//  printf("[%d]END GCMCoupler::call_ice_model() calling model->run_timestep()\n", rank);

    // Record what the ice model produced.
    // NOTE: This shouldn't change model->ice_ovals_I
    IceModel_Writer *owriter = model->owriter();
    if (owriter) {
        printf("BEGIN owriter->run_timestep()\n");
        owriter->run_decoded(time_s, model->ice_ovals_I);
        printf("END owriter->run_timestep()\n");
    }

printf("END GCMCoupler::call_ice_model(nfields=%ld)\n", nfields);
};

/** @param sbuf the (filled) array of ice grid values for this MPI node. */
void GCMCoupler::couple_to_ice(
double time_s,
int nfields,            // Number of fields in sbuf.  Not all will necessarily be filled, in the case of heterogeneous ice models.
ibmisc::DynArray<SMBMsg> &sbuf, // Values, already converted to ice model inputs (from gcm outputs)
std::vector<SparseVector> &gcm_ivals)   // Root node only: Already-allocated space to put output values.  Members as defined by the VarSet GCMCoupler::gcm_inputs
{
    // TODO: Convert this to use ibmisc::gather_msg_array() instead!!!

    // Gather buffers on root node
    int num_mpi_nodes;
    MPI_Comm_size(gcm_params.gcm_comm, &num_mpi_nodes); 

printf("BEGIN GCMCoupler::couple_to_ice() time_s=%f, sbuf.size=%d, sbuf.ele_size=%d\n", time_s, sbuf.size, sbuf.ele_size);

    // MPI_Gather the count
    std::unique_ptr<int[]> rcounts;
    rcounts.reset(new int[num_mpi_nodes]);
    for (int i=0; i<num_mpi_nodes; ++i) rcounts[i] = 0;

    int nele_l = sbuf.size;
    MPI_Gather(&nele_l, 1, MPI_INT, &rcounts[0], 1, MPI_INT, gcm_params.gcm_root, gcm_params.gcm_comm);

    // Compute displacements as prefix sum of rcounts
    std::unique_ptr<int[]> displs;
    std::unique_ptr<ibmisc::DynArray<SMBMsg>> rbuf;

    displs.reset(new int[num_mpi_nodes+1]);
    displs[0] = 0;
    for (int i=0; i<num_mpi_nodes; ++i) displs[i+1] = displs[i] + rcounts[i];
    int nele_g = displs[num_mpi_nodes];

    // Create receive buffer, and gather into it
    // (There's an extra item in the array for a sentinel)
    rbuf.reset(new ibmisc::DynArray<SMBMsg>(SMBMsg::size(nfields), nele_g+1));

    MPI_Datatype mpi_type(SMBMsg::new_MPI_struct(nfields));
    MPI_Gatherv(sbuf.begin().get(), sbuf.size, mpi_type,
        rbuf->begin().get(), &rcounts[0], &displs[0], mpi_type,
        gcm_params.gcm_root, gcm_params.gcm_comm);
    MPI_Type_free(&mpi_type);

    if (am_i_root()) {
        // (ONLY ON GCM ROOT)
        // Clear output arrays, which will be filled in additively
        // on each ice model
//      for (auto ov=gcm_ivals.begin(); ov != gcm_ivals.end(); ++ov) *ov = 0;

        // Add a sentinel
        (*rbuf)[rbuf->size-1].sheetix = 999999;

        // Sort the receive buffer so items in same ice sheet
        // are found together
        qsort(rbuf->begin().get(), rbuf->size, rbuf->ele_size, &SMBMsg::compar);

        // (ONLY ON GCM ROOT)
        // Figure out which ice sheets we have data for
        auto lscan(rbuf->begin());
        auto rscan(lscan);
        std::map<int, CallIceModelParams> im_params;
        while (rscan < rbuf->end()) {
            if (rscan->sheetix != lscan->sheetix) {
                int sheetix = lscan->sheetix;
                auto cimp(CallIceModelParams(sheetix, lscan.get(), rscan.get()));
                im_params[sheetix] = cimp;
                lscan = rscan;
            }

            ++rscan;
        }

        // (ONLY ON GCM ROOT)
        // NOTE: call_ice_model() is called (below) even on NON-ROOT
        // Call all our ice models
        for (size_t sheetix=0; sheetix < ice_models.size(); ++sheetix) {
            auto &ice_model(ice_models[sheetix]);

            // Assume we have data for all ice models
            // (So we can easily maintain MPI SIMD operation)
            auto params(im_params.find(sheetix));
            call_ice_model(&*ice_model, sheetix, time_s, *rbuf,
                params->second.begin, params->second.next);

            // Convert to variables the GCM wants (but still on the ice grid)
            ice_model->set_gcm_inputs(0);   // Fills in gcm_ivals_I

            // Free ice_ovals_I
            ice_model->free_ice_ovals_I();
        }

        regrid_gcm_inputs_onroot(time_s, gcm_ivals, 0);
    } else {
        // (ONLY ON NOT GCM ROOT)
        // We're not root --- we have no data to send to ice
        // models, we just call through anyway because we will
        // receive data in an upcomming MPI_Scatter
        // Call all our ice models
        for (size_t sheetix=0; sheetix < ice_models.size(); ++sheetix) {
            auto &ice_model(ice_models[sheetix]);

            // Assume we have data for all ice models
            // (So we can easily maintain MPI SIMD operation)
            call_ice_model(&*ice_model, sheetix, time_s, *rbuf,
                NULL, NULL);

        }       // if (gcm_params.gcm_rank == gcm_params.gcm_root)

    }

    printf("END GCMCoupler::couple_to_ice()\n");
}


void GCMCoupler::scaled_regrids(
    std::string const regrid_spec,
    std::vector<SparseMatrix> AvIs,
    SparseVector scaleA)
{
    for (int sheetix=0; sheetix < ice_models.size(); ++sheetix) {
        RegridMatrices rm(&*regridder.sheets[sheetix]);
        auto AvI(rm.regrid(regrid_spec, false));

        AvIs.push_back(std::move(AvI->M));

        for (auto ii=AvI->weight.begin(); ii != AvI->weight.end(); ++ii)
            scaleA.add(*ii, ii.val());
    }

    scaleA.consolidate({0});
    for (auto ii=scaleA.begin(); ii != scaleA.end(); ++ii)
        ii.val() = 1. / ii.val();
}


/** Implicit input: gcm_ivals_I, set by set_gcm_inputs() */
void GCMCoupler::regrid_gcm_inputs_onroot(
double time_s,
std::vector<SparseVector> &gcm_ivals,   // Root node only: Already-allocated space to put output values.  Members as defined by the VarSet GCMCoupler::gcm_inputs
unsigned int mask)
{
    // This method is meant to be run only on the GCM root.
    if (!am_i_root()) return;

    printf("BEGIN GCMCoupler::regrid_gcm_inputs_onroot\n");

    // (ONLY ON GCM ROOT)
    // =============== Regrid to the grid requested by the GCM

    // (ONLY ON GCM ROOT)
    // ----------- Get regrid matrices for A<-I and E<-I
    std::vector<SparseMatrix> AvIs;
    SparseVector scaleA;
    scaled_regrids("AvI", AvIs, scaleA);

    std::vector<SparseMatrix> EvIs;
    SparseVector scaleE;
    scaled_regrids("EvI", EvIs, scaleE);


    // (ONLY ON GCM ROOT)
    // ------------ Regrid each GCM input from ice grid to whatever grid it needs.
    for (int var_ix=0; var_ix < gcm_inputs.size(); ++var_ix) {
        VarMeta const &cf(gcm_inputs.data[var_ix]);

        if ((cf.flags & mask) != mask) continue;

        auto dest_grid(cf.flags & contracts::GRID_BITS);

        // --- Assemble all inputs, to multiply by ice_to_hp matrix

        // Accumulate sum of all ice sheets for this variable
        SparseVector ival_unscaled;
        for (size_t sheetix=0; sheetix < ice_models.size(); ++sheetix) {
            auto &ice_model(ice_models[sheetix]);

            // Compute: ival_unscaled += XvI * valI
            auto &XvIs(dest_grid == contracts::ATMOSPHERE ? AvIs : EvIs);
            auto &XvI(XvIs[sheetix]);
            auto &valI(ice_model->gcm_ivals_I[var_ix]);
            for (auto ii=XvI.begin(); ii != XvI.end(); ++ii)
                ival_unscaled.add({ii.index(0)}, ii.val() * valI(ii.index(1)));;
        }
        ival_unscaled.consolidate({0});

        // Scale it!
        auto &scaleX(dest_grid == contracts::ATMOSPHERE ? scaleA : scaleE);
        multiply_ele(gcm_ivals[var_ix], ival_unscaled, scaleX);
    }
    // ----------------- Free Memory
    // (ONLY ON GCM ROOT)
    for (size_t sheetix=0; sheetix < ice_models.size(); ++sheetix) {
        auto &ice_model(ice_models[sheetix]);

        ice_model->free_ovals_ivals_I();
    }

    printf("END GCMCoupler::regrid_gcm_inputs_onroot\n");
}


/** Follows the pattern of couple_to_ice()
@param sbuf the (filled) array of ice grid values for this MPI node. */
void GCMCoupler::get_initial_state(
double time_s,
std::vector<SparseVector> &gcm_ivals)   // Root node only: Already-allocated space to put output values.  Members as defined by the VarSet GCMCoupler::gcm_inputs
{
    printf("BEGIN GCMCoupler::get_initial_state()\n");

    if (am_i_root()) {

        // (ONLY ON GCM ROOT)
        // NOTE: call_ice_model() is called (below) even on NON-ROOT
        // Call all our ice models
        for (size_t sheetix=0; sheetix < ice_models.size(); ++sheetix) {
            auto &ice_model(ice_models[sheetix]);

            ice_model->get_initial_state(time_s);

            // Record what the ice model produced.
            // NOTE: This shouldn't change ice_model->ice_ovals_I
            IceModel_Writer *owriter = ice_model->owriter();
            const double time_s = gcm_params.time_start_s;
            if (owriter) {
                printf("BEGIN owriter->run_timestep()\n");
                owriter->run_decoded(time_s, ice_model->ice_ovals_I);
                printf("END owriter->run_timestep()\n");
            }

            // Convert to variables the GCM wants (but still on the ice grid)
            ice_model->set_gcm_inputs(contracts::INITIAL);  // Fills in gcm_ivals_I

            // Free ice_ovals_I
            ice_model->free_ice_ovals_I();
        }

        // Fill in gcm_ivals
        regrid_gcm_inputs_onroot(time_s, gcm_ivals, contracts::INITIAL);

    } else {
        // (NOT GCM ROOT)
        // We're not root --- we have no data to send to ice
        // models, we just call through anyway because we will
        // receive data in an upcomming MPI_Scatter
        // Call all our ice models
        for (size_t sheetix=0; sheetix < ice_models.size(); ++sheetix) {
            auto &ice_model(ice_models[sheetix]);

            ice_model->get_initial_state(time_s);
        }       // if (gcm_params.gcm_rank == gcm_params.gcm_root)

    }

    printf("END GCMCoupler::get_initial_state()\n");
}

int GCMCoupler::rank() const
    { return gcm_params.gcm_rank; }

}       // namespace
