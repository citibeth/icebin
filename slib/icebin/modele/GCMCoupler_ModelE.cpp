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
#include <ibmisc/string.hpp>
#include <ibmisc/f90blitz.hpp>
#include <icebin/modele/GCMCoupler_ModelE.hpp>
#include <icebin/contracts/contracts.hpp>
#include <icebin/domain_splitter.hpp>
#include <boost/filesystem.hpp>
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


using namespace ibmisc;

namespace icebin {
namespace modele {

static double const nan = std::numeric_limits<double>::quiet_NaN();
using namespace icebin::contracts;

// ======================================================================

/** @param endj Last(+1) element of each domain (0-based indexing) */
DomainDecomposer_ModelE::DomainDecomposer_ModelE(
    std::vector<int> const &endj,
    ibmisc::Domain const &_domainA_global) :    // Starts from ModelE; j indexing base=1
domainA_global(_domainA_global),
rank_of_j(endj[endj.size()-1], blitz::fortranArray)    // Allocate for 1-based indexing
{
    ndomain = endj.size();    // startj contains an extra sentinel item at the end
    int j=0;
    for (int irank=0; irank<ndomain; ++irank) {
        for (; j < endj[irank]; ++j)
            rank_of_j(j) = irank;    // zero-based indexing for j

    }
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

    scalars.add("by_dt", nan, "s-1", 1., "Inverse of coupling timestep");

    scalars.add("unit", nan, "", 0, "Dimensionless identity");
//  gcm_input_scalars.add_field("unit", nan, "", 0, "Dimensionless identity");


}
// -----------------------------------------------------
// Called from LISnow::allocate()

static std::vector<HCSegmentData> parse_hc_segments(std::string const &str)
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


extern "C" void *gcmce_new(
    ModelEParams const &_rdparams,

    // Info about the global grid
    int im, int jm,

    // Info about the local grid (1-based indexing)
    int i0, int i1, int j0, int j1,

    // MPI Stuff
    MPI_Fint comm_f, int root)
{
    std::unique_ptr<GCMCoupler_ModelE> self(
        new GCMCoupler_ModelE());

    ModelEParams &rdparams(self->rdparams);
    rdparams = _rdparams;
    GCMParams &params(self->gcm_params);

    // Domains and indexing are alphabetical indexes, zero-based
    params.domainA = ibmisc::Domain({i0+1,j0+1}, {i1, j1});
    params.domainA_global = ibmisc::Domain({1,1}, {im+1, jm+1});
    params.gcm_comm = MPI_Comm_f2c(comm_f);
    MPI_Comm_rank(params.gcm_comm, &params.gcm_rank);
    self->world.reset(new boost::mpi::communicator(
        params.gcm_comm, boost::mpi::comm_attach));
    params.gcm_root = root;

    params.icebin_config_fname = boost::filesystem::absolute("./ICEBIN_IN").string();
    params.config_dir = boost::filesystem::canonical("./ICEBIN_MODEL_CONFIG_DIR").string();
    params.run_dir = boost::filesystem::absolute(".").string();
    params.gcm_dump_dir = boost::filesystem::absolute("ibdump").string();

    params.hc_segments = parse_hc_segments(f_to_cpp(
        rdparams.icebin_segments, sizeof(rdparams.icebin_segments)));
    params.icebin_base_hc = params.ec_segment().base;



    // Read the coupler, along with ice model proxies
    if (params.am_i_root()) self->ncread(params.icebin_config_fname, "m");

    // Check bounds on the IceSheets, set up any state, etc.
    // This is done AFTER setup of self because self->read_from_netcdf()
    // might change the IceSheet, in certain cases.
    // (for example, if PISM is used, elev2 and mask2 will be read from related
    // PISM input file, and the version in the ICEBIN file will be ignored)


    /** Creates (on root) a picture of the full domain decomposition.
    Run from all MPI ranks... */
    std::vector<int> endj;
    int endme = params.domainA[1].end;
    boost::mpi::gather<int>(*self->world,
        &endme, 1,    // In-values
        endj, 0);                    // Out-values, root
    if (params.am_i_root()) {
        self->domains.reset(new DomainDecomposer_ModelE(endj, params.domainA_global));
    }

    // TODO: Test that im and jm are consistent with the grid read.
    return self.release();
}
// ==========================================================
// Called from LISheetIceBin::read_nhc_gcm()

/* Tells ModelE how many elevation classes it needs **/
extern "C"
int gcmce_nhc_gcm(GCMCoupler_ModelE *self)
    { return self->nhc_gcm(); }

int GCMCoupler_ModelE::read_nhc_gcm()
{
    NcIO ncio(this->gcm_params.icebin_config_fname, 'r');
    int nhc_ice = ncio.nc->getDim("m.nhc").getSize();

    // Find the "ec" segment and set its size now...
    auto &ec(this->gcm_params.ec_segment());
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
void gcmce_add_gcm_outputE(
GCMCoupler_ModelE *self,
F90Array<double, 3> &var_f,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
char const *long_name_f, int long_name_len)
{
    std::string field_name(field_name_f, field_name_len);
    std::string units(units_f, units_len);
    std::string long_name(long_name_f, long_name_len);
    std::unique_ptr<blitz::Array<double,3>> var(
        new blitz::Array<double,3>(var_f.to_blitz()));

    unsigned int flags = 0;

    static double const xnan = std::numeric_limits<double>::quiet_NaN();
    self->gcm_outputsE.add(
        field_name, xnan, units, flags, long_name);

    self->modele_outputs.gcm_ovalsE.push_back(std::move(var));
}
// -----------------------------------------------------
/** @para var_nhp Number of elevation points for this variable.
 (equal to 1 for atmosphere variables, or nhp for elevation-grid variables)
@param return: Start of this variable in the gcm_inputs_local array (Fortran 1-based index) */
extern "C"
void gcmce_add_gcm_inputA(
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
    std::unique_ptr<blitz::Array<double,2>> var(
        new blitz::Array<double,2>(var_f.to_blitz()));

    unsigned int flags = 0;
    if (initial) flags |= contracts::INITIAL;

    static double const xnan = std::numeric_limits<double>::quiet_NaN();
    self->gcm_inputsAE[GridAE::A].add(
        field_name, xnan, units, flags, long_name);

    self->modele_inputs.gcm_ivalsA.push_back(std::move(var));
}
// -----------------------------------------------------
extern "C"
void gcmce_add_gcm_inputE(
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
    std::unique_ptr<blitz::Array<double,3>> var(
        new blitz::Array<double,3>(var_f.to_blitz()));

    unsigned int flags = 0;
    if (initial) flags |= contracts::INITIAL;

    static double const xnan = std::numeric_limits<double>::quiet_NaN();
    self->gcm_inputsAE[GridAE::E].add(
        field_name, xnan, units, flags, long_name);

    self->modele_inputs.gcm_ivalsE.push_back(std::move(var));
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
    auto &modele_inputs(self->modele_inputs);
    modele_inputs.fhc.reference(fhc.to_blitz());
    modele_inputs.elevE.reference(elevE.to_blitz());
    modele_inputs.focean.reference(focean.to_blitz());
    modele_inputs.flake.reference(flake.to_blitz());
    modele_inputs.fgrnd.reference(fgrnd.to_blitz());
    modele_inputs.fgice.reference(fgice.to_blitz());
    modele_inputs.zatmo.reference(zatmo.to_blitz());
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
void gcmce_couple_native(GCMCoupler_ModelE *self,
int itime,
bool run_ice);    // if false, only initialize

extern "C"
void gcmce_cold_start(GCMCoupler_ModelE *self, int itimei)
{
    // This will cold-start initial conditions of the dynamic ice model

    // Call superclass cold_start()
    self->cold_start(
        ibmisc::Datetime(self->rdparams.yeari,1,1),
        itimei * self->rdparams.dtsrc);

    // b) Compute fhc, elevE
    // c) Compute ZATMO, FGICE, etc.
    self->update_topo();
    if (self->gcm_params.dynamic_topo) {
        // TODO...
    }

    // d) Sync with dynamic ice model
    gcmce_couple_native(self, itimei, false);

    // Inits files used to dump gcm_in and gcm_out
    if (self->gcm_params.gcm_dump_dir.size() > 0)
        boost::filesystem::create_directory(self->gcm_params.gcm_dump_dir);
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
bool run_ice)    // if false, only initialize
{
    double time_s = itime * self->rdparams.dtsrc;

    // Fill it in...
    VectorMultivec gcm_ovalsE_s(self->gcm_outputsE.size());

    auto &indexingA(self->gcm_regridder.gridA->indexing);

    // domain uses alphabetical order, 0-based indexing...
    const auto base_hc(self->gcm_params.icebin_base_hc);
    const auto nhc_ice(self->gcm_regridder.nhc(0));
    std::vector<double> val(self->gcm_outputsE.size());    // Temporary

    auto &domainA(self->domainA);
    for (int ihc=base_hc; ihc < base_hc + nhc_ice; ++ihc) {
    for (int j=domainA[1].begin; j < domainA[1].end; ++j) {
    for (int i=domainA[0].begin; i < domainA[0].end; ++i) {
        // i,j are 0-based indexes.
        int i_f = i+1;
        int j_f = j+1;
        int ihc_f = ihc+1;

        if (self->modele_inputs.fhc(i_f,j_f,ihc_f) == 0) continue;    // C2F indexing

        int ihc_ice = ihc - base_hc - 1;   // Use only IceBin HC's, F2C index conversion

        long iE = self->gcm_regridder.indexingE.tuple_to_index(
            make_array(i, j, ihc_ice));

        for (unsigned int ivar=0; ivar<self->gcm_outputsE.size(); ++ivar) {
            val[ivar] = (*self->modele_outputs.gcm_ovalsE[ivar])(i_f,j_f,ihc_f);    // Fortran-order,
        }

        gcm_ovalsE_s.add({iE}, &val[0]);
    }}}


    // Gather it to root
    // boost::mpi::communicator &gcm_world(world);
    // Init our output struct based on number of A and E variables.
    GCMInput out({self->gcm_inputsAE[0].size(), self->gcm_inputsAE[1].size()});
    if (self->gcm_params.am_i_root()) {
        std::vector<VectorMultivec> every_gcm_ovalsE_s;
        boost::mpi::gather(*self->world, gcm_ovalsE_s,
            every_gcm_ovalsE_s, self->gcm_params.gcm_root);

        // Concatenate coupler inputs
        VectorMultivec gcm_ovalsE_s(concatenate(every_gcm_ovalsE_s));

        // Couple on root!
        GCMInput out(
            self->couple(time_s,
                gcm_ovalsE_s, run_ice));

        // Split up the output (and 
        std::vector<GCMInput> every_outs(
            split_by_domain<DomainDecomposer_ModelE>(out,
                *self->domains, *self->domains));

        // Scatter!
        boost::mpi::scatter(*self->world, every_outs, out, self->gcm_params.gcm_root);


    } else {    // ~world.am_i_root()
        // Send our input to root
        boost::mpi::gather(*self->world, gcm_ovalsE_s, self->gcm_params.gcm_root);

        // Let root do the work...
        // TODO... we must call through MPI so ice model can do MPI scatter/gathers..
        self->couple(time_s

        // Receive our output back from root
        boost::mpi::scatter(*self->world, out, self->gcm_params.gcm_root);
    }

    // 1. Copies values back into modele.gcm_ivals
    self->update_gcm_ivals(out);
    // 2. Sets icebin_nhc, 
    // 3. Updates FHC, ZATMO, etc.
    self->update_topo();
}
// =======================================================
/** Called from MPI rank */
void GCMCoupler_ModelE::update_gcm_ivals(GCMInput const &out)
{
    // Update gcm_ivalA variables...
    // Read from here...
    VectorMultivec const &gcm_ivalsA_s(out.gcm_ivalsAE[GridAE::A]);
    // Write to here...
    std::vector<std::unique_ptr<blitz::Array<double,3>>> &gcm_ivalsA(
        gcm_ivalsA);
    auto nvar(out.nvar());    // A and E

    for (size_t i=0; i<gcm_ivalsA_s.size(); ++i) {
        long iA = gcm_ivalsA_s.index[i];
        auto ij(gcm_regridder.indexing(GridAE::A).index_to_tuple<int,2>(iA));    // zero-based, alphabetical order
        int const i_f = ij[0]+1;    // C2F
        int const j_f = ij[1]+1;
        for (unsigned int ivar=0; ivar<nvar[GridAE::A]; ++ivar) {
            (*gcm_ivalsA[ivar])(i_f, j_f) =
                gcm_ivalsA_s.vals[i*nvar[GridAE::A] + ivar];
        }
    }

    // Update gcm_ivalE variables...
    // Read from here...
    VectorMultivec const &gcm_ivalsE_s(out.gcm_ivalsAE[GridAE::E]);
    // Write to here...
    std::vector<std::unique_ptr<blitz::Array<double,3>>> &gcm_ivalsE(
        gcm_ivalsE);

    for (size_t i=0; i<gcm_ivalsE_s.size(); ++i) {
        long iE = gcm_ivalsE_s.index[i];
        auto ijk(gcm_regridder.indexing(GridAE::E).index_to_tuple<int,2>(iE));
        int const i_f = ijk[0]+1;    // C2F
        int const j_f = ijk[1]+1;
        int const ihc_ice = ijk[2];    // zero-based, just EC's known by ice model
        int const ihc_gcm_f = gcm_params.icebin_base_hc + ihc_ice + 1;    // 1-based, G

        for (unsigned int ivar=0; ivar<nvar[GridAE::E]; ++ivar) {
            (*gcm_ivalsE[ivar])(i_f, j_f, ihc_gcm_f) =
                gcm_ivalsE_s.vals[i*nvar[GridAE::E] + ivar];
        }
    }
}




// NetCDF Logging Stuff





}}
