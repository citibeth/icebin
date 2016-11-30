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

    // domain uses Fortran-order, 0-based indexing...
    for (int ihc=icebin_base_hc; ihc < icebin_base_hc + icebin_nhc; ++ihc) {
    for (int j=domainA.base[1]; j != domainA.base[1]+domainA.extent[1]; ++j) {
    for (int i=domainA.base[0]; i != domainA.base[0]+domainA.extent[0]; ++i) {
        // i,j are 0-based indexes.

        if (modele_f.fhc(i+1,j+1,ihc+1) == 0) continue;    // C2F indexing

        long iE = gcm_regridder->indexingHC.tuple_to_index<2>({
            indexingA.tuple_to_index<2>({i, j}),
            ihc-icebin_base_hc-1   // Use only IceBin HC's, F2C index conversion
        });

        for (unsigned int ivar=0; ivar<nvar; ++ivar)
            val[ivar] = gcm_ovalsE(i,j,ihc,ivar);    // Fortran-order

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

        // Log output from GCM
        write_gcm_ovalsE(last_time_s, gcm_ovalsE_s);

        // every_outs has one item per MPI rank
        std::vector<GCMCouplerOutput> every_outs(
            this->couple(time_s,
                to_array(concatenate(every_gcm_ovalsE_s))));

        // Log input back to GCM
        this->write_gcm_coupler_outputs(time_s, every_outs);

        // Scatter!
        boost::mpi::scatter(world, every_outs, out, world.root);


    } else {    // ~world.am_i_root()
        // Send our input to root
        boost::mpi::gather(world, gcm_ovalsE_s, world.root);

        // Let root do the work...

        // Receive our output back from root
        boost::mpi::scatter(world, out, world.root);
    }

    // Update gcm_ival variables...
    for (unsigned iAE=0; iAE<2; ++iAE) {
        VectorSparseParallelVectors &gcm_ivalsX_s(out.gcm_ivals[iAE]);
        std::vector<blit::Array<double,3>> &gcm_ivalsX(modele.gcm_ivals[iAE]);

        std::array<int,2> ahc;
        std::array<int,2> ij;
        for (size_t i=0; i<out.gcm_ivalsX_s.size(); ++i) {
            if (iAE == GCMI::A) {
                long iA = gcm_ivalsX_s.index[i];
                auto ij(indexingA.index_to_tuple(iA));    // zero-based
                int const i_f = ij[0]+1;    // C2F
                int const j_f = ij[1]+1;
                for (unsigned int ivar=0; ivar<out.nvar; ++ivar) {
                    gcm_ivalsX[ivar](i_f, j_f) =
                        out.gcm_ivalsX_s.vals[i*out.nvar + ivar];
                }

            } else {
                long iE = gcm_ivalsX_s.index[i];
                ahc = indexingHC.index_to_tuple(iE);
                ij = indexingA.index_to_tuple(ahc[0]);
                int const i_f = ij[0]+1;    // C2F
                int const j_f = ij[1]+1;
                int const ihc_f = ahc[1]+icebin_base_hc+1;
                gcm_ivalsX(i_f, j_f, ihc_f) =
                    out.gcm_ivalsX_s.vals[i*out.nvar + ivar];
            }
        }
    }

    // Update FHC, TOPO, etc.
    // TODO: Nothing for now
}
// ===============================================================
// NetCDF Logging Stuff





}}
