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
#include <ibmisc/memory.hpp>
#include <icebin/modele/GCMCoupler_ModelE.hpp>
#include <icebin/contracts/contracts.hpp>

use namespace ibmisc;

namespace icebin {
namespace modele {

static double const nan = std::numeric_limits<double>::quiet_NaN();
using namespace icebin::contracts;

// ---------------------------------------------------------------
GCMCoupler_ModelE::GCMCoupler_ModelE() :
    GCMCoupler(GCMCoupler::Type::MODELE)
{

    // ------------ GCM Outputs
    // The GCM must produce the same set of outputs, no matter what
    // ice model is being used
    gcm_outputs.add_field("runo", nan, "kg m-2", ELEVATION,
        "Downward water flux through bottom layer");
    gcm_outputs.add_field("eruno", nan, "J m-2", ELEVATION,
        "Enthalpy of downward water flux through bottom layer");
#ifdef TRACERS_WATER
    gcm_outputs.add_field("trruno")
#endif
    gcm_outputs.add_field("deltah", nan, "J m-2", ELEVATION,
        "Enthalpy change of 'borrowed' layer");
    gcm_outputs.add_field("massxfer", nan, "kg m-2", ELEVATION,
        "Mass of ice being transferred Stieglitz --> Icebin");
    gcm_outputs.add_field("enthxfer", nan, "J m-2", ELEVATION,
        "Enthlpy of ice being transferred Stieglitz --> Icebin");
#ifdef TRACERS_WATER
    gcm_outputs.add_field("trxfer");
#endif
    gcm_outputs.add_field("volxfer", nan, "m^3 m-2", ELEVATION,
        "Volume of ice being transferred Stieglitz --> Icebin");

//    gcm_outputs.add_field("unit", nan, "", 0, "Dimensionless identity");


    // ------------------------- GCM Inputs
    // ModelE sets this, via repeated calls to add_gcm_input_ij()
    // and add_gcm_input_ijhc().  See alloc_landic_com() in LANDICE_COM.f

    // ----------------- Scalars provided by the GCM
    // Scalars are things that can only be computed at the last minute
    // (eg, dt for a particular coupling timestep).  Constants that
    // can be computed at or before contract initialization time can
    // be placed directly into the VarTransformer.

    scalars.add_field("by_dt", nan, "s-1", 1., "Inverse of coupling timestep");

    scalars.add_field("unit", nan, "", 0, "Dimensionless identity");
//  gcm_input_scalars.add_field("unit", nan, "", 0, "Dimensionless identity");

}
// ---------------------------------------------------------------
std::unique_ptr<GCMPerIceSheetParams>
GCMCoupler_ModelE::read_gcm_per_ice_sheet_params(
    ibmisc::NcIO &ncio,
    std::string const &sheet_vname)
{


    // Read GCM-specific coupling parameters
    // Set the contract for each ice sheet, based on:
    //   (a) GCM-specific coupling parameters (to be read),
    //   (b) The type of ice model

    auto gcm_var = giss::get_var_safe(nc, (sheet_vname + ".modele").c_str());

    std::unique_ptr<GCMPerIceSheetParams_ModelE> params(
        new GCMPerIceSheetParams_ModelE());

    params->coupling_type = giss::parse_enum<ModelE_CouplingType>(
        giss::get_att(gcm_var, "coupling_type")->as_string(0));

    return static_cast_unique_ptr<GCMPerIceSheetParams>(params);
}
// ---------------------------------------------------------------
struct ModelEMsg {
    int iA, ihp;    // Indices into ModelE
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
/**
This will:
1. Simultaneously:
   a) Subtract 1 (Fotran base) off of all indices
   b) Subtract hc_offset from ihc index of each non-zero element
   c) Convert indexes to single iE index
   d) MPI Gather
   e) Convert to sparse form on MPI root node.

2. Call couple()

3. Reverse (1)

Put elsewhere in the Fortran code:

3. Update FHC, zatmo, etc. based on new elevE and AvE

4. Update state variables based on new elevation grid (new ice sheet elevations)


@param gcm_ovalsE(nhc_gcm,jm,im) Dense arrays directly from ModelE:
   a) 1-based indexing
   b) Distributed across multiple MPI nodes
   c) ihc is ihc_gcm

@param hc_offset Offset of first IceBin elevation class in GCM's set
    of EC's (zero-based).
*/
void GCMCoupler_ModelE::couple_native(
double time_s,
blitz::Array<double,4> const &gcm_ovalsE,    // Fortran array: im,jm,nhc,nvar
bool do_run)
{
//    int im = gcm_ovalsE.extent(0);
//    int jm = gcm_ovalsE.extent(1);

    int nhc_gcm = gcm_ovalsE.extent(2);
    int nvar = gcm_ovalsE.extent(3);

    // Count total number of elements in the inputs (for this MPI domain)
    long nele_l = 0;
    for (int ihc=icebin_base_hc; ihc < icebin_base_hc + icebin_nhc; ++ihc) {
    for (int j=domainA.base[0]; j != domainA.base[0]+domainA.extent[0]; ++j) {
    for (int i=domainA.base[1]; i != domainA.base[1]+domainA.extent[1]; ++i) {
        if (fhc(im,jm,nhc) != 0) ++nele_l;
    }}}

    // Allocate buffer for that amount of stuff
    ibmisc::DynArray<ModelEMsg> sbuf(ModelEMsg::size(nvar), nele_l);

    // Fill it in...
    std::vector<int> tupleA;
    for (int ihc=icebin_base_hc; ihc < icebin_base_hc + icebin_nhc; ++ihc) {
    for (int j=domainA.base[0]; j != domainA.base[0]+domainA.extent[0]; ++j) {
    for (int i=domainA.base[1]; i != domainA.base[1]+domainA.extent[1]; ++i) {
        if (fhc(im,jm,nhc) == 0) continue;

        ModelEMsg &msg(sbuf[nmsg]);
        tupleA[0] = j-1;    // F2C index conversion
        tupleA[1] = i-1;    // F2C index conversion
        msg.iA = indexingA.tuple_to_index(tupleA);
        msg.ihc = ihc-icebin_base_hc-1;    // Use only IceBin HC's, F2C index conversion
        for (unsigned int ivar=0; ivar<nvar; ++ivar)
            msg[ivar] = gcm_ovalsE(i,j,ihc,ivar);
        ++nmsg;
    }}}

    // Sanity check: make sure we haven't overrun our buffer
    if (nmsg != sbuf.size) (*icebin_error)(
        "Wrong number of items in buffer: %d vs %d expected\n", nmsg, sbuf.size);

    // Gather it to root
    GCMParams const &gcm_params(api->gcm_coupler.gcm_params);
    std::unique_ptr<ibmisc::DynArray<ModelEMsg>> rbuf = ibmisc::gather_msg_array(
        gcm_params.gcm_comm, gcm_params.gcm_root, sbuf, nfields, nmsg, 0);

    // Process the gathered data
    if (rank == gcm_params.gcm_root) {

        SparseParallelVectorsE gcm_ovalsE_s(rbuf->size());

        // Construct ArraySparseParallelVectorsE from msgs
        size_t msg_size = ModelEMsg::size(nvar);    // bytes
        gcm_ovalsE_s.ixA.reference(make_array(
            &rbuf[0].iA, msg_size / sizeof(long)));
        gcm_ovalsE_s.ixHC.reference(make_array(
            &rbuf[0].ihc, msg_size / sizeof(int)));
        gcm_ovalsE_s.values.reserve(nvar);
        for (int ivar=0; ivar<nvar; ++ivar) {
            gcm_ovalsE_s.values.push_back(
                make_array(&rbuf[0].vals[i], msg_size / sizeof(double)));
        }

        // Couple!
        this->couple(time_s, gcm_ovalsE_s, do_run);

    }

    // Get dimensions of full domain
    int nhp_gcm = icebin_modele_nhp_gcm(api);
    ModelEDomain const *domain(&*api->domain);

}


}}
