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
#include <icebin/domain_splitter.hpp>
#include <boost/algorithm/string.hpp>

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

// ======================================================================

/** @param endj Last(+1) element of each domain (0-based indexing) */
DomainDecomposer_ModelE::DomainDecomposer_ModelE(std::vector<int> const &lastj, im_world, jm_world) :    // Starts from ModelE; j indexing base=0
    rank_of_j(endj[endj.size()]),    // 0-based indexing
    im_world(_im_world), jm_world(_jm_world);
{
    ndomain = endj.size();    // startj contains an extra sentinel item at the end
    int j=0;
    for (int irank=0; irank<ndomain; ++irank) {
        for (; j < endj[irank]; ++j)
            rank_of_j(j) = irank;    // zero-based indexing for j
}

/** Creates (on root) a picture of the full domain decomposition.
Called from all MPI ranks... */
std::unique_ptr<DomainDecomposer_ModelE> new_domain_decomposer_mpi(
    ibmisc::Domain const &domainA_global,    // alphabetical order, zero-based indexing
    ibmisc::Domain const &domainA)
{
    std::vector<int> endj;
    boost::mpi::gather(world, domainA.end[1], endj, 0);

    std::unique_ptr<DomainDecomposer_ModelE> ret;
    if (am_i_root) {
        ret.reset(new DomainDecomposer_ModelE(endj, domainA_global.high[0], domainA_global.high[1]));
    }
    return ret;
}


// ======================================================================
// ======================================================================
// ======================================================================
// Called from LISnow::allocate()
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
// -----------------------------------------------------
// Called from LISnow::allocate()
extern "C" void *gcmce_new(
    ModelEParams const &_rdparams,

    // Info about the global grid
    int im, int jm,

    // Info about the local grid (1-based indexing)
    int i0, int i1, int j0, int j1,

    // MPI Stuff
    MPI_Fint comm_f, int root)
{
    std::unique_ptr<GCMCoupler_ModelE> gcm_coupler(
        new GCMCoupler_ModelE());

    ModelEParams &rdparams(gcm_coupler->rdparams);
    rdparams = rdparams;
    GCMParams &params(gcm_coupler->gcm_params);

    // Domains and indexing are alphabetical indexes, zero-based
    params.domainA = ibmisc::Domain({i0+1,j0+1}, {i1, j1});
    params.domainA_global = ibmisc::Domain({1,1}, {im+1, jm+1});
    params.gcm_comm = MPI_Comm_f2c(comm_f);
    gcm_coupler->world.reference(new boost::mpi::communicator(
        gcm_params.comm, boost::mpi::comm_attach));
    params.gcm_root = root;

    params.icebin_config_fname = boost::filesystem::absolute("./ICEBIN_IN");
    params.config_dir = boost::filesystem::canonical("./ICEBIN_MODEL_CONFIG_DIR");
    params.run_dir = boost::filesystem::absolute(".");

    params.hc_segments = parse_hc_segments(f_to_cpp(
        rundeck.icebin_segments, sizeof(rundeck.icebin_segments)));

    // Read the coupler, along with ice model proxies
    gcm_coupler->ncread(gcm_params.icebin_config_fname, "m", gcm_params.domainA);

    // Check bounds on the IceSheets, set up any state, etc.
    // This is done AFTER setup of gcm_coupler because self->read_from_netcdf()
    // might change the IceSheet, in certain cases.
    // (for example, if PISM is used, elev2 and mask2 will be read from related
    // PISM input file, and the version in the ICEBIN file will be ignored)

    gcm_coupler.domains = std::move(new_domain_decomposer_mpi(domainA_global, domainA));

    // TODO: Test that im and jm are consistent with the grid read.
    return gcm_coupler.release();
}

static std::string parse_hc_segments(std::string const &str)
{
    std::vector<HCSegmentData> ret;

    // TDOO: Parse for real later
    std::vector<std::string> segment_names;
    boost::algorithm::split(segment_names, str, boost::is_any_of(","));
    int base = 0;
    for (auto &seg : segment_names) {
        boost::algorithm::trim(seg);
        int len;
        if (seg == "legacy") len = 1;
        else if (seg == "sealand") len = 2;
        else len = -1;

        ret.push_back(HCSegmentData(seg, base, len));
        base += len;
    }
    return ret;
}
// ==========================================================
// Called from LISheetIceBin::read_nhc_gcm()

/* Tells ModelE how many elevation classes it needs **/
extern "C"
int gcmce_read_nhc_gcm(GCMCoupler_ModelE *self)
{
    NcIO ncio(self->gcm_params.icebin_config_fname, 'r');
    int nhc_ice = ncio.nc->getDim("m.nhc").getSize();

    // Find the "ec" segment and set its size now...
    auto &ec(self->gcm_params.ec_segment());
    ec.size = nhc_ice;
    return ec.base + ec.size;

    // NcIO Destructor closes...
}
// ---------------------------------------------------------------
// ==========================================================
// Called from LISheetIceBin::allocate()

/** Set a single constant value in Icebin.  This is a callback, to be called
from ModelE's (Fortran code) constant_set::set_all_constants() */
/** Called from within LISheetIceBin::allocate() */
void gcmce_set_constant(
    GCMCoupler_ModelE *self,
    char const *name_f, int name_len,
    double val,
    char const *units_f, int units_len,
    char const *description_f, int description_len)
{
    self->gcm_constants.set(
        std::string(name_f, name_len),
        val,
        std::string(units_f, units_len),
        std::string(description_f, description_len));
}

// ---------------------------------------------------------------
extern "C"
int gcmce_add_gcm_outputE(
GCMCoupler_ModelE *self,
F90Array<double, 3> &var_f,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
char const *long_name_f, int long_name_len)
{
    std::string field_name(field_name_f, field_name_len);
    std::string units(units_f, units_len);
    std::string long_name(long_name_f, long_name_len);
    auto var(var_f.to_blitz());

    unsigned int flags = 0;

    static double const xnan = std::numeric_limits<double>::quiet_NaN();
    self->gcm_outputs.add(
        field_name, xnan, units, flags, long_name);

    self->modele_outputs.gcm_ovalsE.push_back(var);

    printf("gcmce_add_gcm_inputA(%s, %s, %s) --> %d\n", field_name.c_str(), units.c_str(), grid.c_str(), ret);

    return ret;
}
// -----------------------------------------------------
/** @para var_nhp Number of elevation points for this variable.
 (equal to 1 for atmosphere variables, or nhp for elevation-grid variables)
@param return: Start of this variable in the gcm_inputs_local array (Fortran 1-based index) */
extern "C"
int gcmce_add_gcm_inputA(
GCMCoupler_ModelE *self,
F90Array<double, 2> &var_f,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
int initial,    // bool
char const *long_name_f, int long_name_len)
{
    std::string field_name(field_name_f, field_name_len);
    std::string units(units_f, units_len);
    std::string long_name(long_name_f, long_name_len);
    auto var(new_unique_ptr(var_f.to_blitz());

    unsigned int flags = 0;
    if (initial) flags |= contracts::INITIAL;

    static double const xnan = std::numeric_limits<double>::quiet_NaN();
    self->gcm_inputs[GridAE::A].add(
        field_name, xnan, units, flags, long_name);

    self->modele_inputs.gcm_ivals[GridAE::A].push_back(var);

    printf("gcmce_add_gcm_inputA(%s, %s, %s) --> %d\n", field_name.c_str(), units.c_str(), ret);

    return ret;
}
// -----------------------------------------------------
extern "C"
int gcmce_add_gcm_inputE(
GCMCoupler_ModelE *self,
F90Array<double, 3> &var_f,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
int initial,    // bool
char const *long_name_f, int long_name_len)
{
    std::string field_name(field_name_f, field_name_len);
    std::string units(units_f, units_len);
    std::string long_name(long_name_f, long_name_len);
    auto var(var_f.to_blitz());

    unsigned int flags = 0;
    if (initial) flags |= contracts::INITIAL;

    static double const xnan = std::numeric_limits<double>::quiet_NaN();
    self->gcm_inputs[GridAE::E].add(
        field_name, xnan, units, flags, long_name);

    self->modele_inputs.gcm_ivals[GridAE::E].push_back(var);

    printf("gcmce_add_gcm_inputE(%s, %s, %s) --> %d\n", field_name.c_str(), units.c_str(), ret);

    return ret;
}
// -----------------------------------------------------
// ==========================================================
// Called from LIShetIceBin::reference_globals()

extern "C"
void gcmce_reference_globals(
    GCMCoupler_ModelE *self,
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

// ===========================================================
// Called from LISheetIceBin::io_rsf()   (warm starts)

/** Only called from MPI Root */
extern "C"
void gcmce_io_rsf(GCMCoupler_ModelE *self,
    char *fname_c, int fname_n)
{
    std::string fname(fname_c, fname_n);
    // TODO: Get ice model to save/restore state
    // We must figure out how to get an appropriate filename
}



// ===========================================================
// Called from LISheetIceBin::cold_start()

extern "C"
void gcmce_cold_start(GCMCoupler_ModelE *self, int itimei, double dtsrc)
{
    // a) Cold-start initial conditions of dynamic ice model
TODO...

    // b) Compute fhc, elevE
    // c) Compute ZATMO, FGICE, etc.
    if (self->gcm_params.dynamic_topo) {
        // TODO...
    }

    // d) Sync with dynamic ice model
    gcmce_couple_native(self, itime, dtsrc, false);

    // Inits files used to dump gcm_in and gcm_out
    if (self->gcm_params.gcm_dump_dir.size() > 0)
        mkdir(self->gcm_params.gcm_dump_dir);
}

// =======================================================
// Called from LISheetIceBin::couple()
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
extern "C"
void gcmce_couple_native(GCMCoupler_ModelE *self,
int itime,
int yy, int mm, int dd, // Date that time_s lies on
//std::array<int,3> const &yymmdd, // Date that time_s lies on
bool run_ice)    // if false, only initialize
{
    double time_s = itime * dtsrc;

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

        // Concatenate coupler inputs
        ArraySparseParallelVectors gcm_ovalsE_s(
            to_array(concatenate(every_gcm_ovalsE_s)));

        // Couple on root!
        GCMCouplerOutput out(
            this->couple(time_s, yymmdd,
                gcm_ovalsE_s, run_ice));

        // Split up the output (and 
        std::vector<GCMCouplerOutput<3>> every_outs(split_by_domain(out, domains));

        // Scatter!
        boost::mpi::scatter(world, every_outs, out, world.root);


    } else {    // ~world.am_i_root()
        // Send our input to root
        boost::mpi::gather(world, gcm_ovalsE_s, world.root);

        // Let root do the work...

        // Receive our output back from root
        boost::mpi::scatter(world, out, world.root);
    }

    // 1. Copies values back into modele.gcm_ivals
    modele.update_gcm_ivals(out);
    // 2. Sets icebin_nhc, 
    // 3. Updates FHC, ZATMO, etc.
    modele.update_modele_vars(out);
}
// =======================================================
// =======================================================
// =======================================================

#if 0
/** Callback from new_ice_coupler() */
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
#endif
// -----------------------------------------------------
// ===============================================================

/** Called from MPI rank */
void ModelEInputs::update_gcm_ivals(GCMCouplerOutput const &out)
{
    // Update ModelE variables (zatmo, fhc, etc);
    // also icebin_base_hc and icebin_nhc
    update_modele_vars(out);

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
}




// NetCDF Logging Stuff





}}
