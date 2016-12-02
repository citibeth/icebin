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

// ==========================================================
/** @param nc The IceBin configuration file */
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
    auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});

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

    ice_couplers.clear();
    for (size_t i=0; i < regridder.sheets.size(); ++i) {
        IceRegridder *sheet = &*regridder.sheets[i];
        std::string vname_sheet = vname + "." + sheet->name();

        // Create an IceCoupler corresponding to this IceSheet.
        std::unique_ptr<IceCoupler> ice_coupler(new_ice_coupler(ncio, vname_sheet, this, sheet));
        ice_coupler->ncread(ncio, vname_sheet);

        contracts::setup(*this, *ice_coupler);    // where does this go w.r.t ncread() and upate?
        ice_coupler->update_ice_sheet(ncio, vname_sheet);

        // ---- Create the affiliated writers
        // Writer for ice model input
        auto iwriter(dynamic_cast_unique_ptr<IceCoupler_Writer, IceCoupler>(
            new_ice_coupler(IceCoupler::Type::WRITER, this, sheet)));
        iwriter->init(IceCoupler::INPUT, &*ice_coupler);        // Writer-specific initialization
        ice_coupler->_iwriter = std::move(iwriter);

        // Writer for ice model output
        auto owriter(dynamic_cast_unique_ptr<IceCoupler_Writer, IceCoupler>(
            new_ice_coupler(IceCoupler::Type::WRITER, this, sheet)));
        owriter->init(IceCoupler::OUTPUT, &*ice_coupler);       // Writer-specific initialization
        ice_coupler->_owriter = std::move(owriter);

        ice_couplers.push_back(std::move(ice_coupler));
    }

    ncio.close();
    printf("END GCMCoupler::read_from_netcdf()\n");
}



void GCMCoupler::set_start_time(
    ibmisc::time::tm const &time_base,
    double time_start_s)
{
    gcm_params.set_start_time(time_base, time_start_s);

    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        auto &ice_coupler(ice_couplers[sheetix]);

        ice_coupler->start_time_set();

        // Handle writing stuff, if the
        IceCoupler_Writer *iwriter = ice_coupler->iwriter();    // The affiliated input-writer (if it exists).
        if (iwriter) iwriter->start_time_set();

        IceCoupler_Writer *owriter = ice_coupler->owriter();
        if (owriter) owriter->start_time_set();

        // This function is the last phase of initialization.
        // Only now can we assume that the contracts are fully set up.
#if 1
        // Print out the contract and var transformations
        std::cout << "========= Contract for " << ice_coupler->name() << std::endl;
        std::cout << "---- GCM->Ice     Output Variables:" << std::endl;
        std::cout << ice_coupler->contract[IceCoupler::INPUT];
        std::cout << "TRANSFORMATIONS:" << std::endl;
        std::cout << ice_coupler->var_transformer[IceCoupler::INPUT];
        std::cout << "---- Ice->GCM     Output Variables:" << std::endl;
        std::cout << ice_coupler->contract[IceCoupler::OUTPUT];
        std::cout << "TRANSFORMATIONS:" << std::endl;
        std::cout << ice_coupler->var_transformer[IceCoupler::OUTPUT];
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

/** PROTECTED method.  Called by all MPI nodes. */
void GCMCoupler::call_ice_coupler(
    IceCoupler *model,
    int sheetix,
    double time_s,
    ibmisc::DynArray<SMBMsg> &rbuf,
    SMBMsg *begin, SMBMsg *end)     // Messages have the MAX number of fields for any ice model contract
{
    // ----------------- Construct input arrays
    // The number of fields for THIS ice sheet will be <= the number
    // of fields in SMBMsg
    int nfields = model->contract[IceCoupler::INPUT].size();

printf("BEGIN GCMCoupler::call_ice_coupler(%s, nfields=%ld)\n", model->name().c_str(), nfields);

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
    IceCoupler_Writer *iwriter = model->iwriter();
    if (iwriter) iwriter->run_timestep(time_s, indices, ivals2);

    // Now call to the ice model
//  printf("[%d] BEGIN GCMCoupler::call_ice_coupler() calling model->run_timestep()\n", rank);
    model->run_timestep(time_s, indices, ivals2);
//  printf("[%d]END GCMCoupler::call_ice_coupler() calling model->run_timestep()\n", rank);

    // Record what the ice model produced.
    // NOTE: This shouldn't change model->ice_ovals_I
    IceCoupler_Writer *owriter = model->owriter();
    if (owriter) {
        printf("BEGIN owriter->run_timestep()\n");
        owriter->run_decoded(time_s, model->ice_ovals_I);
        printf("END owriter->run_timestep()\n");
    }

printf("END GCMCoupler::call_ice_coupler(nfields=%ld)\n", nfields);
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
        std::map<int, CallIceCouplerParams> im_params;
        while (rscan < rbuf->end()) {
            if (rscan->sheetix != lscan->sheetix) {
                int sheetix = lscan->sheetix;
                auto cimp(CallIceCouplerParams(sheetix, lscan.get(), rscan.get()));
                im_params[sheetix] = cimp;
                lscan = rscan;
            }

            ++rscan;
        }

        // (ONLY ON GCM ROOT)
        // NOTE: call_ice_coupler() is called (below) even on NON-ROOT
        // Call all our ice models
        for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
            auto &ice_coupler(ice_couplers[sheetix]);

            // Assume we have data for all ice models
            // (So we can easily maintain MPI SIMD operation)
            auto params(im_params.find(sheetix));
            call_ice_coupler(&*ice_coupler, sheetix, time_s, *rbuf,
                params->second.begin, params->second.next);

            // Convert to variables the GCM wants (but still on the ice grid)
            ice_coupler->set_gcm_inputs(0);   // Fills in gcm_ivals_I

            // Free ice_ovals_I
            ice_coupler->free_ice_ovals_I();
        }

        regrid_gcm_inputs_onroot(time_s, gcm_ivals, 0);
    } else {
        // (ONLY ON NOT GCM ROOT)
        // We're not root --- we have no data to send to ice
        // models, we just call through anyway because we will
        // receive data in an upcomming MPI_Scatter
        // Call all our ice models
        for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
            auto &ice_coupler(ice_couplers[sheetix]);

            // Assume we have data for all ice models
            // (So we can easily maintain MPI SIMD operation)
            call_ice_coupler(&*ice_coupler, sheetix, time_s, *rbuf,
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
    for (int sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        RegridMatrices rm(&*regridder.sheets[sheetix]);
        auto AvI(rm.regrid(regrid_spec, false, true));

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
        for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
            auto &ice_coupler(ice_couplers[sheetix]);

            // Compute: ival_unscaled += XvI * valI
            auto &XvIs(dest_grid == contracts::ATMOSPHERE ? AvIs : EvIs);
            auto &XvI(XvIs[sheetix]);
            auto &valI(ice_coupler->gcm_ivals_I[var_ix]);
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
    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        auto &ice_coupler(ice_couplers[sheetix]);

        ice_coupler->free_ovals_ivals_I();
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
        // NOTE: call_ice_coupler() is called (below) even on NON-ROOT
        // Call all our ice models
        for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
            auto &ice_coupler(ice_couplers[sheetix]);

            ice_coupler->get_initial_state(time_s);

            // Record what the ice model produced.
            // NOTE: This shouldn't change ice_coupler->ice_ovals_I
            IceCoupler_Writer *owriter = ice_coupler->owriter();
            const double time_s = gcm_params.time_start_s;
            if (owriter) {
                printf("BEGIN owriter->run_timestep()\n");
                owriter->run_decoded(time_s, ice_coupler->ice_ovals_I);
                printf("END owriter->run_timestep()\n");
            }

            // Convert to variables the GCM wants (but still on the ice grid)
            ice_coupler->set_gcm_inputs(contracts::INITIAL);  // Fills in gcm_ivals_I

            // Free ice_ovals_I
            ice_coupler->free_ice_ovals_I();
        }

        // Fill in gcm_ivals
        regrid_gcm_inputs_onroot(time_s, gcm_ivals, contracts::INITIAL);

    } else {
        // (NOT GCM ROOT)
        // We're not root --- we have no data to send to ice
        // models, we just call through anyway because we will
        // receive data in an upcomming MPI_Scatter
        // Call all our ice models
        for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
            auto &ice_coupler(ice_couplers[sheetix]);

            ice_coupler->get_initial_state(time_s);
        }       // if (gcm_params.gcm_rank == gcm_params.gcm_root)

    }

    printf("END GCMCoupler::get_initial_state()\n");
}

int GCMCoupler::rank() const
    { return gcm_params.gcm_rank; }


// =============================================================

/** Updates regridding matrices using the latest state of the ice
    models. */
void GCMCoupler::update_regridders()
{
    for (auto &ii : ice_regridders) {
        ii->
    }
}


struct SparseParallelVectors {
    // Stores a bunch of parallel sparse vectors
    // Index of each element in the parallel vectors
    std::vector<long> index;
    // Values of for each element in the vectors.  
    std::vector<double> vals;
    // Number of _vals element per _ix element
    int nvar;
};

struct GCMCoupleOutput {
    // Outputs from IceCoupler, transformed and regridded back to E/A
//    std::vector<SparseVector> gcm_ivals;    // both A and E grids.


    // Mapping from the index of a variable in gcm_ivalsE/gcm_ivalsA
    // and the index within the GCMCoupler::gcm_inputs
    SparseParallelVectors[GCMCoupler::GCMI::COUNT] gcm_ivals; // gcm_ivalsE, gcm_ivalsA


    // Values required to update TOPO, etc. in ModelE
    // (see add_fhc.py for how these are to be used)
    // We can get these from AvE
    // SparseVector wAvE;     // Area of A grid cells that overlap ice
    //SparseVector areaA;    // Total (native) area of A grid cells
    //SparseVector elevA;    // Used for atmosphere orography (ZATMO in ModelE)

    // Regrid matrix to go from last step's elevation classes to this
    // step's elevation classes.
    SparseMatrix E1vE0;

    // Regrid matrix to convert to atmosphere.
    // (EvA is assumed by GCM, as long as AvE is local; see Fischer&Nowicki 2014)
    WeightedSparse AvE;

    // Used for temperature downscaling according to a lapse rate
    SparseVector elevE;
};

// ------------------------------------------------------------
// ------------------------------------------------------------

CoupleReturn GCMCoupler::couple(
// Simulation time [s]
double time_s,
// Values from GCM, passed GCM -> Ice
std::vector<SparseVector> gcm_ovalsE,
CoupleReturn &ret)
{
    
    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        auto &ice_coupler(ice_couplers[sheetix]);

        // Store regridding matrices for the last timestep, which we will
        // need to create ice_ivals
        std::unique_ptr<WeightedSparse> IvE0(std::move(ice_coupler.IvE));
        std::unique_ptr<WeightedSparse> IvA0(std::move(ice_coupler.IvA));

        // Create ice_ivals (dense)
        auto ice_ivalsI(allocate_dense(ice_coupler.contract[

        std::vector<blitz::Array<double,1>> ice_ivalsI(gcm_inputs.size());

        VarSet const &ocontract(contract[IceCoupler::OUTPUT]);
        for (int i=0; i < ocontract.size(); ++i)
            ice_ovals_I.push_back(blitz::Array<double,1>(ndata()));



    for (int i=0; i < gcm_inputs.size(); ++i) {
        gcm_ivals_I.push_back(blitz::Array<double,1>(ndata()));
    }



        // Update elevI and compute IvE for later
        RegridMatrices rm(ice_coupler.get_elevI());

        WeightedSparse const &IvE(ice_coupler.IvE);
        WeightedSparse &IvA(rm.regrd("IvA"));




        // Get the matrices we need to create ice_ivals
        WeightedSparse EvI(rm.regrid("EvI", true, true));
        WeightedSparse AvI(rm.regrid("AvI", true, true));


        multiply(ret.Ev1Ev0, EvI.M, IvE0.M);


    }
    return ret;
}


}       // namespace
