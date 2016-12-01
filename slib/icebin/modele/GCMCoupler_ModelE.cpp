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
#include <ibmisc/ncfile.hpp>
#include <icebin/modele/GCMCoupler_ModelE.hpp>
#include <icebin/contracts/contracts.hpp>

// See here to serialize objects with non-default constructor
//    http://www.boost.org/doc/libs/1_62_0/libs/serialization/doc/serialization.html#constructors
// http://www.boost.org/doc/libs/1_62_0/libs/serialization/doc/tutorial.html
// http://www.ocoudert.com/blog/2011/07/09/a-practical-guide-to-c-serialization/
// Binary archive that defines boost::archive::binary_oarchive
// and boost::archive::binary_iarchive
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>


use namespace ibmisc;

namespace icebin {
namespace modele {

static double const nan = std::numeric_limits<double>::quiet_NaN();
using namespace icebin::contracts;

// ---------------------------------------------------------------
GCMCoupler_ModelE::GCMCoupler_ModelE() :
    GCMCoupler(GCMCoupler::Type::MODELE)
{

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
std::array<int,3> const &yymmdd, // Date that time_s lies on
bool init_only)
{
    GCMCouplerOutput out;
//    int im = gcm_ovalsE.extent(0);
//    int jm = gcm_ovalsE.extent(1);

    int nhc_gcm = gcm_ovalsE.extent(2);
    int nvar = gcm_ovalsE.extent(3);

    // Allocate buffer for that amount of stuff
    ibmisc::DynArray<ModelEMsg> sbuf(ModelEMsg::size(nvar), nele_l);


    // Fill it in...
    VectorSparseParallelVectors gcm_ovalsE_s(nvar);
    std::vector<double> val(nvar);

    auto &indexingA(regridder.gridA->indexing);

    // domain uses alphabetical order, 0-based indexing...
    for (int ihc=icebin_base_hc; ihc < icebin_base_hc + icebin_nhc; ++ihc) {
    for (int j=domainA.base[1]; j != domainA.base[1]+domainA.extent[1]; ++j) {
    for (int i=domainA.base[0]; i != domainA.base[0]+domainA.extent[0]; ++i) {
        // i,j are 0-based indexes.
        int i_f = i+1
        int j_f = j+1
        int ihc_f = ihc+1

        if (modele_f.fhc(i_f,j_f,ihc_f) == 0) continue;    // C2F indexing

        int ihc_ice = ihc-icebin_base_hc-1   // Use only IceBin HC's, F2C index conversion
        long iE = gcm_regridder->indexingHC.tuple_to_index<2>({
            indexingA.tuple_to_index<2>({i, j}),
            ihc_ice
        });

        for (unsigned int ivar=0; ivar<nvar; ++ivar) {
            val[ivar] = (*modele_output.gcm_ovalsE[ivar])(i_f,j_f,ihc_f);    // Fortran-order,
        }

        gcm_ovalsE_s.add(iE, val);
    }}}


    // Gather it to root
    // boost::mpi::communicator &gcm_world(world);
    GCMCouplerOutput out;
    if (world.am_i_root()) {
        std::vector<VectorSparseParallelVectors> every_gcm_ovalsE_s;
        boost::mpi::gather(world, gcm_ovalsE_s, every_gcm_ovalsE_s, world.root);

        // Concatenate them all
        ArraySparseParallelVectors gcm_ovalsE_s(
            to_array(concatenate(every_gcm_ovalsE_s)));

        // every_outs has one item per MPI rank
        GCMCouplerOutput out(
            this->couple(time_s, yymmdd,
                gcm_ovalsE_s, init_only));

        // Scatter!
        boost::mpi::scatter(world, every_outs, out, world.root);


    } else {    // ~world.am_i_root()
        // Send our input to root
        boost::mpi::gather(world, gcm_ovalsE_s, world.root);

        // Let root do the work...

        // Receive our output back from root
        boost::mpi::scatter(world, out, world.root);
    }

    // Update gcm_ivalA variables...
    VectorSparseParallelVectors &gcm_ivalsA_s(out.gcm_ivals[GridAE::A]);
    std::vector<std::unique_ptr<blitz::Array<double,3>>> &gcm_ivalsA(modele.gcm_ivals[GridAE::A]);

    for (size_t i=0; i<out.gcm_ivalsA_s.size(); ++i) {
        long iA = gcm_ivalsA_s.index[i];
        auto ij(indexingA.index_to_tuple(iA));    // zero-based, alphabetical order
        int const i_f = ij[0]+1;    // C2F
        int const j_f = ij[1]+1;
        for (unsigned int ivar=0; ivar<out.nvar; ++ivar) {
            (*gcm_ivalsA[ivar])(i_f, j_f) =
                out.gcm_ivalsA_s.vals[i*out.nvar + ivar];
        }
    }

    // Update gcm_ivalE variables...
    VectorSparseParallelVectors &gcm_ivalsE_s(out.gcm_ivals[GridAE::E]);
    std::vector<std::unique_ptr<blitz::Array<double,3>>> &gcm_ivalsE(modele.gcm_ivals[GridAE::E]);

    for (size_t i=0; i<out.gcm_ivalsE_s.size(); ++i) {
        long iE = gcm_ivalsE_s.index[i];
        ijk = indexingE.index_to_tuple(iE);
        int const i_f = ijk[0]+1;    // C2F
        int const j_f = ijk[1]+1;
        int const ihc_ice = ijk[2];    // zero-based, just EC's known by ice model
        int const ihc_gcm_f = icebin_base + ihc_ice + 1;    // 1-based, G

        for (unsigned int ivar=0; ivar<out.nvar; ++ivar) {
            (*gcm_ivalsE[ivar])(i_f, j_f, ihc_gcm_f) =
                out.gcm_ivalsE_s.vals[i*out.nvar + ivar];
        }
    }

    // Update FHC, TOPO, etc.
    // TODO: Nothing for now
}
// ===============================================================
// NetCDF Logging Stuff





}}
