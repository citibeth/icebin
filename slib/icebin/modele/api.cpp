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
#include <ibmisc/netcdf.hpp>
//#include <ibmisc/mpi.hpp>
#include <ibmisc/blitz.hpp>
#include <ibmisc/f90blitz.hpp>
//#include <icebin/HCIndex.hpp>
#include <icebin/modele/icebin_modele.hpp>
#include <icebin/modele/GCMCoupler_ModelE.hpp>
//#include <icebin/IceModel_TConv.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <icebin/contracts/contracts.hpp>
#include <ibmisc/exit.hpp>

using namespace icebin;
using namespace icebin::modele;
using namespace ibmisc;
using namespace netCDF;


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


// ================================================================

struct ModelEMsg {
    int i, j, k;    // Indices into ModelE
    double vals[1];     // Always at least one val; but this could be extended, based on # of inputs

    double &operator[](int i) { return *(vals + i); }

    /** @return size of the struct, given a certain number of values */
    static size_t size(int nfields)
        { return sizeof(ModelEMsg) + (nfields-1) * sizeof(double); }

    static MPI_Datatype new_MPI_struct(int nfields);

    /** for use with qsort */
//  static int compar(void const * a, void const * b);

};

MPI_Datatype ModelEMsg::new_MPI_struct(int nfields)
{
    int nele = 3 + nfields;
    int blocklengths[] = {1, 1, 1, nfields};
    MPI_Aint displacements[] = {offsetof(ModelEMsg,i), offsetof(ModelEMsg,j), offsetof(ModelEMsg,k), offsetof(ModelEMsg, vals)};
    MPI_Datatype types[] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE};
    MPI_Datatype ret;
    MPI_Type_create_struct(4, blocklengths, displacements, types, &ret);
    MPI_Type_commit(&ret);
    return ret;
}
// -----------------------------------------------------

/** LOCAL FUNCTION: Set up a netCDF file, ready to store timeseries
variables in by Icebin.  This is used for both gcm_input and gcm_ouput
files. */
extern "C"
void init_ncfile(icebin_modele *api,
ibmisc::CouplingContract const &contract,
std::string const &fname)
{
    printf("BEGIN init_ncfile(%s)\n", fname.c_str());
    GCMParams const &gcm_params(api->gcm_coupler.gcm_params);

    // Set up NetCDF file to store GCM output as we received them (modele_out.nc)
    int nhp_gcm = icebin_modele_nhp_gcm(api);
    NcIO ncio(fname, NcFile::replace)

    auto dims_elevation(get_or_add_dims(ncio, 
        {"time", "nhp", "jm", "im"}, 
        {-1, nhp_gcm, api->domain->jm, api->domain->im}));
    auto dims_atomsphere(get_or_add_dims(ncio,
        {"time", "jm", "im"}, 
        {-1, api->domain->jm, api->domain->im}));
    auto dims_b(get_or_add_dims(ncio, {"one"}, {1}));

    NcVar grid_var = ncio.nc->addVar("grid", ibmisc::get_nc_type<double>(), dims_b);
    grid_var.putAtt("file", api->gcm_coupler.fname);
    grid_var.putAtt("variable", api->gcm_coupler.vname);

    NcVar time0_var = ncio.nc->addVar("time0", ibmisc::get_nc_type<double>(), dims_b);
    time0_var.putAtt("units", gcm_params.time_units);
    time0_var.putAtt("calendar", "365_day");
    time0_var.putAtt("axis", "T");
    time0_var.putAtt("long_name", "Simulation start time");

    NcVar time_var = ncio.nc->addVar("time", ibmisc::get_nc_type<double>(), dims_atmosphere);
    time_var.putAtt("units", gcm_params.time_units);
    time_var.putAtt("calendar", "365_day");
    time_var.putAtt("axis", "T");
    time_var.putAtt("long_name", "Coupling times");

    //ibmisc::CouplingContract &gcm_outputs(api->gcm_coupler.gcm_outputs);

    for (unsigned int i=0; i < contract.size_nounit(); ++i) {

        NcVar nc_var;
        ibmisc::CoupledField const &cf(contract.field(i));
printf("Creating NC variable for %s (%s)\n", cf.name.c_str(), contracts::to_str(cf.flags).c_str());
        switch(cf.flags & contracts::GRID_BITS) {
            case contracts::ATMOSPHERE :
                nc_var = ncio.nc->addVar(cf.name,
                    ibmisc::get_nc_type<double>(), dims_atmosphere);
            break;
            case contracts::ELEVATION :
                nc_var = ncio.nc->addVar(cf.name,
                ibmisc::get_nc_type<double>(), dims_elevation);
            break;
            default:
                fprintf(stderr, "init_ncfile() unable to handle grid type %s for field %s\n", contracts::to_str(cf.flags), cf.name);
                ibmisc::exit(1);
        }

        auto comment(boost::format(
            "%s[t,...] holds the mean from time[t-1] to time[t].  See time0[0] if t=0.")
            % contract.name(i));
        nc_var.putAtt("comment", comment.str());

        std::string const &description(contract.field(i).get_description());
        if (description != "") nc_var.putAtt("long_name", description);

        std::string const &units(contract.field(i).get_units());
        if (units != "") nc_var.putAtt("units", units);
    }

    // Put initial time in it...
    time0_var.putVar({0}, {1}, &gcm_params.time_start_s);

    ncout.close();
    printf("END init_ncfile(%s)\n", fname.c_str());

}
// -----------------------------------------------------

/** LOCAL FUNCTION: Save the ice model output, after it's been
transformed by Icebin (now it's GCM input). */
static void save_gcm_inputs(
icebin_modele *api,
double time_s,
ibmisc::F90Array<double,3> &gcm_inputs_d_f) // Fortran array
{
    printf("BEGIN save_gcm_inputs(%s)\n", api->gcm_coupler.gcm_in_file.c_str());

    // Fortran-style array: i,j,ihp, indexing starts at 1
    blitz::Array<double,3> gcm_inputs_d(gcm_inputs_d_f.to_blitz());

    // Get dimensions of full domain
    int nhp_gcm = icebin_modele_nhp_gcm(api);
    GCMCoupler &coupler(api->gcm_coupler);
    ibmisc::CouplingContract const &contract(coupler.gcm_inputs);

    // Open output netCDF file
    NcFile ncout(api->gcm_coupler.gcm_in_file.c_str(), NcFile::Write);  // Read/Write
    NcDim *time_dim = ncout.get_dim("time");
    NcDim *nhp_dim = ncout.get_dim("nhp");
    NcDim *jm_dim = ncout.get_dim("jm");
    NcDim *im_dim = ncout.get_dim("im");

    long cur_ijhc[4]{time_dim->size(),0,0,0};       // time, nhp_gcm, jm, im
    long counts_ijhc[4]{1, nhp_dim->size(), jm_dim->size(), im_dim->size()};

    long cur_ij[4]{time_dim->size(),0,0};       // time, nhp, jm, im
    long counts_ij[4]{1, jm_dim->size(), im_dim->size()};

    NcVar *time_var = ncout.get_var("time");
    time_var->set_cur(cur_ijhc);    // Could be cur_ij, either way is fine
    time_var->put(&time_s, counts_ijhc);

    // Write arrays to it
    int base_index = 0;
    for (unsigned int i=0; i < contract.size_nounit(); ++i) {
        double const *array_base = &gcm_inputs_d(1,1,base_index);   // i,j,ihp

        NcVar *nc_var = ncout.get_var(contract.name(i).c_str());

        switch(contract.field(i).flags & contracts::GRID_BITS) {
            case contracts::ATMOSPHERE :
                nc_var->set_cur(cur_ij);
                nc_var->put(array_base, counts_ij);
                base_index += 1;
            break;
            case contracts::ELEVATION :
                nc_var->set_cur(cur_ijhc);
                nc_var->put(array_base, counts_ijhc);
                base_index += nhp_gcm;
            break;
            default: ;
        }
    }

    ncout.close();
    printf("END save_gcm_inputs(%s)\n", api->gcm_coupler.gcm_in_file.c_str());
}
// -----------------------------------------------------

/** LOCAL function: Log the GCM output, which is to be passed to the
ice model.
@param hpvals Values on height-points GCM grid for various fields
    the GCM has decided to provide. */
static void save_gcm_outputs(
icebin_modele *api,
double time_s,
std::vector<std::unique_ptr<blitz::Array<double,3>>> &inputs)
{

    // Get dimensions of full domain
    int nhp_gcm = icebin_modele_nhp_gcm(api);
    ModelEDomain const *domain(&*api->domain);

    int const rank = api->gcm_coupler.rank();   // MPI rank; debugging

    printf("[%d] BEGIN save_gcm_outputs(time_s=%f)\n", rank, time_s);

    GCMCoupler &coupler(api->gcm_coupler);

    ibmisc::CouplingContract &gcm_outputs(coupler.gcm_outputs);

    // Count total number of elements in the inputs (for this MPI domain)
    auto &input0(inputs[0]);
    int nele_l =
        (input0->ubound(2) - input0->lbound(2) + 1) *
        (domain->j1_f - domain->j0_f + 1) *
        (domain->i1_f - domain->i0_f + 1);

    // Find the max. number of fields (for input) used for any ice sheet.
    // This will determine the size of our MPI messages.
    int nfields = gcm_outputs.size_nounit();

    // Allocate buffer for that amount of stuff
    ibmisc::DynArray<ModelEMsg> sbuf(ModelEMsg::size(nfields), nele_l);

    // Fill it in....
    int nmsg = 0;
    for (int k=input0->lbound(2); k<=input0->ubound(2); ++k)        // nhp_gcm
    for (int j=domain->j0_f; j <= domain->j1_f; ++j)
    for (int i=domain->i0_f; i <= domain->i1_f; ++i) {
        ModelEMsg &msg = sbuf[nmsg];
        msg.i = i;
        msg.j = j;
        msg.k = k;
        for (unsigned int l=0; l<nfields; ++l) msg[l] = (*inputs[l])(i,j,k);
        ++nmsg;
    }

    // Sanity check: make sure we haven't overrun our buffer
    if (nmsg != sbuf.size) {
        fprintf(stderr, "Wrong number of items in buffer: %d vs %d expected\n", nmsg, sbuf.size);
        ibmisc::exit(1);
    }

    // Gather it to root
    GCMParams const &gcm_params(api->gcm_coupler.gcm_params);
    std::unique_ptr<ibmisc::DynArray<ModelEMsg>> rbuf = ibmisc::gather_msg_array(
        gcm_params.gcm_comm, gcm_params.gcm_root, sbuf, nfields, nmsg, 0);

    // Process the gathered data
    if (rank == gcm_params.gcm_root) {

        // Allocate ijk arrays
        std::vector<blitz::Array<double,3>> outputs;
        for (unsigned int i=0; i<nfields; ++i) {
            outputs.push_back(blitz::Array<double,3>(nhp_gcm, domain->jm, domain->im));
        }

        // Turn messages into ijk arrays
        for (auto msg=rbuf->begin(); msg != rbuf->end(); ++msg) {
            for (unsigned int i=0; i<nfields; ++i) {
                int mi = msg->i - 1;
                int mj = msg->j - 1;
                int mk = msg->k - 1;        // Convert Fortran --> C indexing
                outputs[i](mk, mj, mi) = (*msg)[i];
            }
        }

        // Write the arrays to a file
        NcFile ncout(api->gcm_coupler.gcm_out_file.c_str(), NcFile::Write); // Read/Write
        NcDim *time_dim = ncout.get_dim("time");
        NcDim *nhp_dim = ncout.get_dim("nhp");
        NcDim *jm_dim = ncout.get_dim("jm");
        NcDim *im_dim = ncout.get_dim("im");

        long cur[4]{time_dim->size(),0,0,0};        // time, nhp_gcm, jm, im
        long counts[4]{1, nhp_dim->size(), jm_dim->size(), im_dim->size()};

        NcVar *time_var = ncout.get_var("time");
        time_var->set_cur(cur);
        time_var->put(&time_s, counts);

        for (int i=0; i<nfields; ++i) {
            NcVar *nc_var = ncout.get_var(gcm_outputs.name(i).c_str());
            nc_var->set_cur(cur);
            nc_var->put(outputs[i].data(), counts);
        }

        ncout.close();
    }
    printf("[%d] END save_gcm_outputs(time_s=%f)\n", rank, time_s);
}

// ================================================================
extern "C" icebin_modele *new_icebin_modele_c()
{
    std::unique_ptr<icebin_modele> api(new icebin_modele());

    // No exception was thrown... we can release our pointer back to Fortran
    icebin_modele *ret = api.release();

//  int const rank = ret->gcm_coupler.rank();   // MPI rank; debugging
//  printf("[%d] Allocated icebin_modele api struct: %p\n", rank, ret);
    return ret;
}
// ---------------------------------------------------
// -----------------------------------------------------
extern "C" void icebin_modele_delete(icebin_modele *&api)
{
    if (api) delete api;
    api = 0;
}
// -----------------------------------------------------
extern "C"
void icebin_modele_reference_modele_vars(
    F90Array<double, 3> &fhc,
    F90Array<double, 3> &elevE,
    F90Array<double, 2> &focean,
    F90Array<double, 2> &flake,
    F90Array<double, 2> &fgrnd,
    F90Array<double, 2> &fgice,
    F90Array<double, 2> &zatmo)
{
    auto &modele(api->gcm_coupler.modele);
    modele.fhc.reference(fhc.to_blitz());
    modele.elevI.reference(elevI.to_blitz());
    modele.focean.reference(focean.to_blitz());
    modele.flake.reference(flake.to_blitz());
    modele.fgrnd.reference(fgrnd.to_blitz());
    modele.fgice.reference(fgice.to_blitz());
    modele.zatmo.reference(zatmo.to_blitz());
}
// -----------------------------------------------------
extern "C"
void icebin_modele_set_start_time(icebin_modele *api, int iyear1, int itimei, double dtsrc)
{
    std::cout << "========= GCM Inputs (second time: must be set by now)" << std::endl;
    std::cout << api->gcm_coupler.gcm_inputs << std::endl;

    GCMParams &gcm_params(api->gcm_coupler.gcm_params);

    api->dtsrc = dtsrc;
    double time0_s = itimei * api->dtsrc;

    api->gcm_coupler.set_start_time(
        ibmisc::time::tm(iyear1, 1, 1),
        time0_s);

    api->itime_last = itimei;

    // Finish initialization...
    // -------------------------------------------------
    // Open file to record GCM outputs to Icebin
    if (api->gcm_coupler.gcm_out_file.length() > 0) {
        init_ncfile(api,
            api->gcm_coupler.gcm_outputs,
            api->gcm_coupler.gcm_out_file);
    }

    // Open file to record GCM inputs from Icebin
    if (api->gcm_coupler.gcm_in_file.length() > 0) {
        init_ncfile(api,
            api->gcm_coupler.gcm_inputs,
            api->gcm_coupler.gcm_in_file);
    }
}
// -----------------------------------------------------
extern "C"
void icebin_modele_get_flice_im_c(icebin_modele *api,
    ibmisc::F90Array<double, 2> &flice1_im_f)       // OUT
{

printf("BEGIN icebin_modele_get_flice_im_c()\n");
std::cout << "flice1_im_f: " << flice1_im_f << std::endl;

    ModelEDomain &domain(*api->domain);

    // Reconstruct arrays, using Fortran conventions
    // (smallest stride first, whatever-based indexing it came with)

    // Get the sparse vector values
    ibmisc::VectorSparseVector<int,double> flice1_s;
    api->gcm_coupler.maker->fgice(flice1_s);

    // Translate the sparse vectors to the ModelE data structures
    std::vector<std::tuple<int, int, double>> flice1_vals;
    for (auto ii = flice1_s.begin(); ii != flice1_s.end(); ++ii) {
        int i1 = ii->first;

        // Filter out things not in our domain
        // (we'll get the answer for our halo via a halo update)
        // Convert to local (ModelE 2-D) indexing convention
        int lindex[domain.num_local_indices];
        domain.global_to_local(i1, lindex);
        if (!domain.in_domain(lindex)) continue;

        // Store it away
        // (we've eliminated duplicates, so += isn't needed, but doesn't hurt either)
        flice1_vals.push_back(std::make_tuple(lindex[0], lindex[1], ii->second));
    }

    // Zero out the ICEBIN-only version completely
    auto flice1_im(flice1_im_f.to_blitz());
    flice1_im = 0;

    // Replace with our values
    for (auto ii=flice1_vals.begin(); ii != flice1_vals.end(); ++ii) {
        int ix_i = std::get<0>(*ii);
        int ix_j = std::get<1>(*ii);
        double val = std::get<2>(*ii);

        flice1_im(ix_i, ix_j) += val;
    }
    // -----------------------------------------------------
    // Balance flice against other landcover types
//  auto fgrnd1(fgrnd1_f.to_blitz());
//  auto focean1(focean1_f.to_blitz());
//  auto flake1(flake1_f.to_blitz());
// This correction is taken care of in ModelE (for now)
// See: FLUXES.f alloc_fluxes()
//  fgrnd1 = 1.0 - focean1 - flake1 - flice1;

#if 0
NcFile nc("flice1_1.nc", NcFile::Replace);
auto flice1_c(ibmisc::f_to_c(flice1));
auto flice1_im_c(ibmisc::f_to_c(flice1_im));
auto a(ibmisc::netcdf_define(nc, "flice1", flice1_c));
auto b(ibmisc::netcdf_define(nc, "flice1_im", flice1_im_c));
a();
b();
nc.close();
#endif

    // -----------------------------------------------------
printf("END icebin_modele_get_flice_im_c()\n");
}
// -----------------------------------------------------
/** Produces the (dense) FHC_IM array from the (sparse) hp_to_atm
coming from raw IceBin. */
extern "C"
void icebin_modele_get_fhc_im_c(icebin::modele::icebin_modele *api,
    ibmisc::F90Array<double, 3> &fhc_im1h_f)    // OUT
{
    auto fhc_im1h(fhc_im1h_f.to_blitz());

    auto &regridder(api->gcm_coupler.regridder);
    auto &indexingHP(regridder.indexingHP);
    // HCIndex &hc_index(*api->gcm_coupler.maker->hc_index);

    RegridMatrices rm(regridder.$$$$$$$$$$$$$$$$$$$$$$
    std::unique_ptr<ibmisc::VectorSparseMatrix> hp_to_atm(regridder.hp_to_atm());
    ModelEDomain &domain(*api->domain);

    // Filter this array, and convert to fhc format
    for (auto ii = hp_to_atm->begin(); ii != hp_to_atm->end(); ++ii) {
        int lindex1a[domain.num_local_indices];

        // Input: HP space
        int lindex[domain.num_local_indices];
        int hp1b, i1b;
        int i3b = ii.col();
        hc_index.index_to_ik(i3b, i1b, hp1b);
        domain.global_to_local(i1b, lindex);
        if (!domain.in_domain(lindex)) {
            //printf("Not in domain: i3b=%d (%d, %d, %d)\n", i3b, lindex[0], lindex[1], hp1b);
            continue;
        }

        // Output: GCM grid
        int i1a = ii.row();
        if (i1a != i1b) {
            fprintf(stderr, "HP2ATM matrix is non-local!\n");
            ibmisc::exit(1);
        }

        // Now fill in FHC_IM
        // +1 for C-to-Fortran conversion
        fhc_im1h(lindex[0], lindex[1], hp1b+1) += ii.val();
    }
    hp_to_atm.release();


    // In a perfect world, summing FHC over elevation points will
    // sum to one.  But in reality, it sums to something else, depending
    // on size difference between grids on sphere vs. on the plane.

}
// -----------------------------------------------------
static void densify_gcm_inputs_onroot(icebin_modele *api,
    std::vector<ibmisc::VectorSparseVector<int,double>> const &gcm_ivals_global,
    ibmisc::F90Array<double,3> &gcm_inputs_d_f)
{
    // This should only be run on the MPI root node
    if (!api->gcm_coupler.am_i_root()) return;

    printf("BEGIN densify_gcm_inputs_onroot()\n");

    int const rank = api->gcm_coupler.rank();   // MPI rank; debugging
    ibmisc::CouplingContract const &contract(api->gcm_coupler.gcm_inputs);
    int n1 = api->gcm_coupler.maker->n1();

    // Fortran-style array: i,j,ihp, indexing starts at 1
    blitz::Array<double,3> gcm_inputs_d(gcm_inputs_d_f.to_blitz());

    // We ARE the root note --- densify the data into the global gcm_inputs array
    for (long ix = 0; ix < contract.size_nounit(); ++ix) {
        int ihp = api->gcm_inputs_ihp[ix];              // First elevation point for this variable
        int var_nhp = api->gcm_inputs_ihp[ix+1] - ihp;  // # elevation points for this variable (in ModelE)

        // Check bounds
        if (ihp+var_nhp > gcm_inputs_d.extent(2)) {
            fprintf(stderr, "[%d] gcm_inputs_d[nhp=%d] is too small (needs at least %d)\n", rank, gcm_inputs_d.extent(2), api->gcm_inputs_ihp[contract.size_nounit()]); //ihp+var_nhp);
            ibmisc::exit(1);
        }

        // Ignore elevation point = 0 (for ELEVATION grid only),
        // which is reserved for ModelE's "legacy" elevation point.
        int icebin_ihp;         // First elevation point for this variable in Icebin
        int icebin_var_nhp; // Number of elevation points for this variable in Icebin
        if (var_nhp == 1) {
            // ATMOSPHERE grid destination, only one "elevation point"
            icebin_ihp = ihp;
            icebin_var_nhp = var_nhp;
        } else {
            // We have an ELEVATION grid destination (not ATMOSPHERE)

            // Skip over elevation point zero for further regridding
            icebin_ihp = ihp + 1;
            icebin_var_nhp = var_nhp - 1;

            // Clear elevation point 0, since it's not involved in the
            // Icebin computation.
            blitz::Array<double,1> ep_zero(
                &gcm_inputs_d(1,1,ihp), // i,j,ihp
                blitz::shape(n1*1), blitz::neverDeleteData);
//          ep_zero = contract.field(ix).default_value;
            ep_zero = 0;
        }

        // Index into our big array-of-array of all gcm_inputs
        // (ignoring the modelE elevation point 0)
        blitz::Array<double,1> dense1d(
            &gcm_inputs_d(1,1,icebin_ihp),  // i,j,ihp
            blitz::shape(n1*icebin_var_nhp), blitz::neverDeleteData);

        // Convert this sparse vector...
        printf("Setting gcm_input %s to %g\n", contract.field(ix).name.c_str(), contract.field(ix).default_value);
//      dense1d = contract.field(ix).default_value;
        dense1d = 0;
        for (auto ii=gcm_ivals_global[ix].begin(); ii != gcm_ivals_global[ix].end(); ++ii) {
            int const i1 = ii->first;
            double const val = ii->second;
            dense1d(i1) += val;
        }
    }

    printf("END densify_gcm_inputs_onroot()\n");


}
// -------------------------------------------------------------
/** @param hpvals Values on height-points GCM grid for various fields
    the GCM has decided to provide.

    @param gcm_inputs_d_f Global (gathered) array of the Icebin
    outputs to be fed into the GCM.  This only needs to be allocated
    if api->gcm_coupler.am_i_root().

    Inputs variables implied by repeated previous calls to icebin_modele_set_gcm_output().
*/
extern "C"
void  icebin_modele_couple_to_ice_c(
icebin_modele *api,
int itime,
ibmisc::F90Array<double,3> &gcm_inputs_d_f)
{
    int rank = api->gcm_coupler.rank(); // MPI rank; debugging
    double time_s = itime * api->dtsrc;

    printf("[%d] BEGIN icebin_modele_couple_to_ice_c(itime=%d, time_s=%f, dtsrc=%f)\n", rank, itime, time_s, api->dtsrc);

    if (api->gcm_coupler.gcm_out_file.length() > 0) {
        // Write out to DESM file
        save_gcm_outputs(api, time_s, api->gcm_outputs);
    }

    // Count total number of elements in the matrices
    // (_l = local to this MPI node)
    int nele_l = 0; //api->gcm_coupler.maker->ice_matrices_size();
printf("icebin_modele_couple_to_ice_c(): hp_to_ices.size() %ld\n", api->hp_to_ices.size());
    for (auto ii = api->hp_to_ices.begin(); ii != api->hp_to_ices.end(); ++ii) {
        nele_l += ii->second.size();
    }

    // Find the max. number of fields (for input) used for any ice sheet.
    // This will determine the size of our MPI messages.
    int nfields_max = 0;
    for (auto sheet=api->gcm_coupler.maker->sheets.begin(); sheet != api->gcm_coupler.maker->sheets.end(); ++sheet) {
        int sheetno = sheet->index;
        ibmisc::VarTransformer &vt(api->gcm_coupler.models[sheetno]->var_transformer[IceModel::INPUT]);
        nfields_max = std::max(nfields_max, (int)vt.dimension(ibmisc::VarTransformer::OUTPUTS).size_nounit());
    }

    // Allocate buffer for that amount of stuff
printf("icebin_modele_couple_to_ice_c(): nfields_max=%d, nele_l = %d\n", nfields_max, nele_l);
    ibmisc::DynArray<SMBMsg> sbuf(SMBMsg::size(nfields_max), nele_l);

    // Fill it in by doing a sparse multiply...
    // (while translating indices to local coordinates)
    HCIndex &hc_index(*api->gcm_coupler.maker->hc_index);
    int nmsg = 0;
printf("[%d] hp_to_ices.size() = %ld\n", rank, api->hp_to_ices.size());
    for (auto ii = api->hp_to_ices.begin(); ii != api->hp_to_ices.end(); ++ii) {
        int sheetno = ii->first;
        ibmisc::VarTransformer &vt(api->gcm_coupler.models[sheetno]->var_transformer[IceModel::INPUT]);


        std::vector<hp_to_ice_rec> &mat(ii->second);

printf("[%d] mat[sheetno=%d].size() == %ld\n", rank, sheetno, mat.size());
        // Skip if we have nothing to do for this ice sheet
        if (mat.size() == 0) continue;

        // Get the CSR sparse matrix to convert GCM outputs to ice model inputs
        ibmisc::CSRAndUnits trans = vt.apply_scalars({
            std::make_pair("by_dt", 1.0 / ((itime - api->itime_last) * api->dtsrc)),
            std::make_pair("unit", 1.0)});

        // Do the multiplication
        for (int j=0; j < mat.size(); ++j) {
            hp_to_ice_rec &jj(mat[j]);
            SMBMsg &msg = sbuf[nmsg];
            msg.sheetno = sheetno;
            msg.i2 = jj.row;

#if 1
            // Convert from (GCM output) to (Ice Model input) units while regridding
            for (int xi=0; xi<vt.dimension(ibmisc::VarTransformer::OUTPUTS).size_nounit(); ++xi) {
                double inp = 0;
                std::vector<std::pair<int, double>> const &row(trans.mat[xi]);
                for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
                    int xj = xjj->first;
                    double io_val = xjj->second;
                    inp += io_val * (*api->gcm_outputs[xj])(jj.col_i, jj.col_j, jj.col_k);
                }
                msg[xi] = jj.val * (inp + trans.units[xi]);
            }
#else
            // This is the old code for the above (that doesn't convert units)
            msg[0] = jj.val * (*api->gcm_outputs[0])(jj.col_i, jj.col_j, jj.col_k);
            msg[1] = jj.val * (*api->gcm_outputs[1])(jj.col_i, jj.col_j, jj.col_k);
            msg[2] = jj.val * (*api->gcm_outputs[2])(jj.col_i, jj.col_j, jj.col_k);
#endif
//          // This is even older code (variables are hardwired)
//          msg[0] = jj.val * smb1h(jj.col_i, jj.col_j, jj.col_k);
//          msg[1] = jj.val * seb1h(jj.col_i, jj.col_j, jj.col_k);
//          msg[2] = jj.val * tg21h(jj.col_i, jj.col_j, jj.col_k);

//printf("msg = %d (i,j, hc)=(%d %d %d) i2=%d %g %g (%g %g)\n", msg.sheetno, lindex[0], lindex[1], ihp+1, msg.i2, msg[0], msg[1], smb1h(lindex[0], lindex[1], ihp+1), seb1h(lindex[0], lindex[1], ihp+1));

            ++nmsg;
        }
    }

    // Sanity check: make sure we haven't overrun our buffer
    if (nmsg != sbuf.size) {
        fprintf(stderr, "Wrong number of items in buffer: %d vs %d expected\n", nmsg, sbuf.size);
        ibmisc::exit(1);
    }

    // sbuf has elements for ALL ice sheets here
    ibmisc::CouplingContract const &contract(api->gcm_coupler.gcm_inputs);
    std::vector<ibmisc::VectorSparseVector<int,double>> gcm_ivals_global(contract.size_nounit());
    api->gcm_coupler.couple_to_ice(time_s, nfields_max, sbuf, gcm_ivals_global);

    // Decode the outputs to a dense array for ModelE
    if (api->gcm_coupler.am_i_root()) {
        // auto gcm_inputs_d(gcm_inputs_d_f.to_blitz());

        densify_gcm_inputs_onroot(api, gcm_ivals_global, gcm_inputs_d_f);

        if (api->gcm_coupler.gcm_in_file.length() > 0) {
            // Write out to DESM file
            save_gcm_inputs(api, time_s, gcm_inputs_d_f);
        }
    }

    api->itime_last = itime;

    printf("[%d] END icebin_modele_couple_to_ice_c(itime=%d)\n", rank, itime);
}
// -------------------------------------------------------------
extern "C"
void  icebin_modele_get_initial_state_c(
icebin_modele *api,
int itime,
ibmisc::F90Array<double,3> &gcm_inputs_d_f)
{

    double time_s = itime * api->dtsrc;

//  auto gcm_inputs_d(gcm_inputs_d_f.to_blitz());

    printf("BEGIN icebin_modele_get_initial_state_c()\n");

    // sbuf has elements for ALL ice sheets here
    ibmisc::CouplingContract const &contract(api->gcm_coupler.gcm_inputs);
    std::vector<ibmisc::VectorSparseVector<int,double>> gcm_ivals_global(contract.size_nounit());
    api->gcm_coupler.get_initial_state(time_s, gcm_ivals_global);

    // Decode the outputs to a dense array for ModelE
    if (api->gcm_coupler.am_i_root()) {

        densify_gcm_inputs_onroot(api, gcm_ivals_global, gcm_inputs_d_f);

        if (api->gcm_coupler.gcm_in_file.length() > 0) {
            save_gcm_inputs(api, time_s, gcm_inputs_d_f);
        }
    }
    api->itime_last = itime;

    printf("END icebin_modele_get_initial_state_c()\n");
}

extern"C"
int icebin_modele_set_gcm_output_c(
icebin_modele *api,
char const *field_name_f, int field_name_len,
ibmisc::F90Array<double, 3> &arr_f)
{
    std::string const field_name(field_name_f, field_name_len);
    int field_ix = api->gcm_coupler.gcm_outputs.index(field_name);
    if (api->gcm_outputs.size() <= field_ix) api->gcm_outputs.resize(field_ix+1);
    blitz::Array<double,3> *arrp = new blitz::Array<double,3>(arr_f.to_blitz());
    api->gcm_outputs[field_ix].reset(arrp);
}

extern "C"
void icebin_modele_couple_native(
    icebin_modele *api,
    double time_s,
    int year, int month, int day,
    bool init_only)
{
    api->gcm_coupler.couple_native(time_s, {year, month, day}, init_only);
}

// ===============================================================
