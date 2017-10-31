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

#include <cstdlib>
#include <mpi.h>        // Intel MPI wants to be first
#include <ibmisc/memory.hpp>
#include <ibmisc/ncfile.hpp>
#include <ibmisc/string.hpp>
#include <ibmisc/f90blitz.hpp>
#include <icebin/modele/GCMCoupler_ModelE.hpp>
#include <icebin/contracts/contracts.hpp>
#include <icebin/domain_splitter.hpp>
#include <boost/filesystem.hpp>
#include <icebin/modele/GCMRegridder_ModelE.hpp>

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

using namespace std;
using namespace ibmisc;
using namespace netCDF;

// See LIGrid.F90
const int UI_ICEBIN = 1;
const int UI_NOTHING= 2;


#if 0
std::ostream &operator<<(std::ostream &os, icebin::modele::ModelEParams const &params)
{
    os << "ModelEParams:" << std::endl;
    os << "    segments=" << std::string(params.icebin_segments,icebin::modele::MAX_CHAR_LEN) << std::endl;
    os << "    dtsrc=" << params.dtsrc << std::endl;
    return os;
}
#endif

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
GCMCoupler_ModelE::GCMCoupler_ModelE(GCMParams &&_params) :
    GCMCoupler(GCMCoupler::Type::MODELE, std::move(_params))
{

    // ----------------- Scalars provided by the GCM
    // Scalars are things that can only be computed at the last minute
    // (eg, dt for a particular coupling timestep).  Constants that
    // can be computed at or before contract initialization time can
    // be placed directly into the VarTransformer.

    scalars.add("by_dt", nan, "s-1", 1., "Inverse of coupling timestep");
}
// -----------------------------------------------------
void GCMCoupler_ModelE::ncread(
    ibmisc::NcIO &ncio_config,
    std::string const &vname)        // comes from this->gcm_params
{
    GCMCoupler::ncread(ncio_config, vname);

    auto config_info(get_or_add_var(ncio_config, vname + ".info", "int64", {}));
    get_or_put_att(config_info, ncio_config.rw, "topo_ocean", topoO_fname);

    // Replace the GCMRegridder with a wrapped version that understands
    // the ocean-vs-atmosphere grid complexity of ModelE
    gcm_regridder.reset(new GCMRegridder_ModelE(gcm_regridder));
}
// -----------------------------------------------------
// Called from LISnow::allocate()

std::string GCMCoupler_ModelE::locate_input_file(
    std::string const &sheet_name,        // eg: greenland
    std::string const &file_name)         // eg: pism_Greenland_5km_v1.1.nc
{
    // Filenames in icebin.cdl have already been absolutized by
    // Modele-Control
    return file_name;
}
// -----------------------------------------------------
// ===========================================================
extern "C"
void *gcmce_new(
    ModelEParams const &_rdparams,

    // Info about the global grid
    int im, int jm,

    // Info about the local grid (1-based indexing, closed ranges)
    int i0_f, int i1_f, int j0_f, int j1_f,

    // MPI Stuff
    MPI_Fint comm_f, int root)
{
printf("BEGIN gcmce_new()\n");

    std::unique_ptr<GCMCoupler_ModelE> self(
        new GCMCoupler_ModelE(GCMParams(MPI_Comm_f2c(comm_f), root)));

    ModelEParams &rdparams(self->rdparams);
    rdparams = _rdparams;
    GCMParams &gcm_params(self->gcm_params);

    // Domains and indexing are alphabetical indexes, zero-based, open ranges
    self->domainA = ibmisc::Domain({i0_f-1,j0_f-1}, {i1_f, j1_f});
    self->domainA_global = ibmisc::Domain({0,0}, {im, jm});

    gcm_params.icebin_config_fname = boost::filesystem::absolute("config/icebin.nc").string();

    // Read the coupler, along with ice model proxies
    self->ncread(gcm_params.icebin_config_fname, "m");

    // Check bounds on the IceSheets, set up any state, etc.
    // This is done AFTER setup of self because self->read_from_netcdf()
    // might change the IceSheet, in certain cases.
    // (for example, if PISM is used, elev2 and mask2 will be read from related
    // PISM input file, and the version in the ICEBIN file will be ignored)


    /** Creates (on root) a picture of the full domain decomposition.
    Run from all MPI ranks... */
    std::vector<int> endj;
    int endme = self->domainA[1].end;
    boost::mpi::gather<int>(self->gcm_params.world,
        &endme, 1,    // In-values
        endj, self->gcm_params.gcm_root);                    // Out-values, root
    if (self->am_i_root()) {
        self->domains.reset(
            new DomainDecomposer_ModelE(endj, self->domainA_global));
    }

    // TODO: Test that im and jm are consistent with the grid read.
    GCMCoupler_ModelE *ret = self.release();
    return ret;
}
// ==========================================================
// Called from LISheetIceBin::_read_nhc_gcm()

/* Tells ModelE how many elevation classes it needs **/
extern "C"
void gcmce_hc_params(GCMCoupler_ModelE *self, int &nhc_gcm, int &icebin_base_hc, int &nhc_ice)
{
    nhc_gcm = self->nhc_gcm();
    icebin_base_hc = self->gcm_params.segment("ec").base;
    nhc_ice = self->gcm_regridder->nhc();
}

int GCMCoupler_ModelE::_read_nhc_gcm()
{
    // Get the name of the grid file
    std::string grid_fname;
    {
        NcIO ncio_config(gcm_params.icebin_config_fname, NcFile::read);
        auto config_info(get_or_add_var(ncio_config, "m.info", "int64", {}));
        get_or_put_att(config_info, ncio_config.rw, "grid", grid_fname);
    }

    // Open the grid file to get the NHC
    NcIO ncio(grid_fname, 'r');
    int nhc_ice = ncio.nc->getDim("m.nhc").getSize();

    // Find the "ec" segment and set its size now...
    auto &ec(this->gcm_params.segment("ec"));
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
extern "C"
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
void gcmce_add_gcm_outpute(
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
/** @para var_nhc Number of elevation points for this variable.
 (equal to 1 for atmosphere variables, or nhc for elevation-grid variables)
@param return: Start of this variable in the gcm_inputs_local array (Fortran 1-based index) */
extern "C"
void gcmce_add_gcm_inputa(
GCMCoupler_ModelE *self,
F90Array<double, 2> &var_f,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
bool initial,    // bool
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
void gcmce_add_gcm_inpute(
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
    F90Array<int, 3> &underice,
    F90Array<double, 3> &elevE,
    F90Array<double, 2> &focean,
    F90Array<double, 2> &flake,
    F90Array<double, 2> &fgrnd,
    F90Array<double, 2> &fgice,
    F90Array<double, 2> &zatmo)
{
    Topos *topos(&self->modele_inputs);
    topos->fhc.reference(fhc.to_blitz());
    topos->underice.reference(underice.to_blitz());
    topos->elevE.reference(elevE.to_blitz());
    topos->focean.reference(focean.to_blitz());
    topos->flake.reference(flake.to_blitz());
    topos->fgrnd.reference(fgrnd.to_blitz());
    topos->fgice.reference(fgice.to_blitz());
    topos->zatmo.reference(zatmo.to_blitz());
}

// ===========================================================
// Called from LISheetIceBin::io_rsf()   (warm starts)

extern "C"
void gcmce_io_rsf(GCMCoupler_ModelE *self,
    char *fname_c, int fname_n)
{
    std::string fname(fname_c, fname_n);
    // TODO: Get ice model to save/restore state
    // We must figure out how to get an appropriate filename
    if (self->am_i_root()) {
    }
}



// ===========================================================
// Called from LISheetIceBin::cold_start()

extern "C"
void gcmce_couple_native(GCMCoupler_ModelE *self,
int itime,
bool run_ice);    // if false, only initialize

extern "C"
void gcmce_cold_start(GCMCoupler_ModelE *self, int yeari, int itimei, double dtsrc)
{
    printf("BEGIN gcmce_cold_start() yeari=%d, itimei=%d, dtsrc=%g\n", yeari, itimei, dtsrc);

    // This will cold-start initial conditions of the dynamic ice model
    self->dtsrc = dtsrc;

    // Call superclass cold_start()
    self->cold_start(
        ibmisc::Datetime(yeari,1,1) ,
        itimei * dtsrc);

    // b) Compute fhc, elevE
    // c) Compute ZATMO, FGICE, etc.
    self->update_topo();
    if (self->gcm_params.dynamic_topo) {
        // TODO...
    }

    // d) Sync with dynamic ice model
    gcmce_couple_native(self, itimei, false);

    printf("END gcmce_cold_start()\n");
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
    double time_s = itime * self->dtsrc;

    // Fill it in...
    VectorMultivec gcm_ovalsE_s(self->gcm_outputsE.size());
    std::vector<double> val(self->gcm_outputsE.size());    // Temporary

    auto &indexingA(self->gcm_regridder->gridA->indexing);

    // domain uses alphabetical order, 0-based indexing...
    const auto base_hc(self->gcm_params.icebin_base_hc);
    const auto nhc_ice(self->gcm_regridder->nhc());

    auto &domainA(self->domainA);
printf("domainA size=%ld base_hc=%d  nhc_ice=%d\n", domainA.data.size(), base_hc, nhc_ice);

    for (int ihc=base_hc; ihc < base_hc + nhc_ice; ++ihc) {
        const int ihc_f = ihc+1;
        const int ihc_ice = ihc - base_hc;   // Use only IceBin HC's, zero-based indexing
        if (ihc_ice < 0) (*icebin_error)(-1,
            "ihc_ice cannot be <0: %d = %d - %d - 1\n", ihc_ice, ihc, base_hc);

        for (int j=domainA[1].begin; j < domainA[1].end; ++j) {
        for (int i=domainA[0].begin; i < domainA[0].end; ++i) {
            // i,j are 0-based indexes.
            const int i_f = i+1;
            const int j_f = j+1;

            if (self->modele_inputs.underice(i_f,j_f,ihc_f) != UI_ICEBIN) continue;    // C2F indexing
            long iE_s = self->gcm_regridder->indexingE.tuple_to_index(
                make_array(i, j, ihc_ice));

            if (iE_s < 0) (*ibmisc_error)(-1,
                "iE_s=%ld (from %d %d %d), it should not be negative\n", iE_s, i, j, ihc_ice);

            for (unsigned int ivar=0; ivar<self->gcm_outputsE.size(); ++ivar) {
                val[ivar] = (*self->modele_outputs.gcm_ovalsE[ivar])(i_f,j_f,ihc_f);    // Fortran-order,
            }

            gcm_ovalsE_s.add({iE_s}, &val[0]);
        }}
    }


    // Gather it to root
    // boost::mpi::communicator &gcm_world(world);
    // Init our output struct based on number of A and E variables.
    GCMInput out({self->gcm_inputsAE[0].size(), self->gcm_inputsAE[1].size()});
    if (self->am_i_root()) {
        // =================== MPI ROOT =============================
        std::vector<VectorMultivec> every_gcm_ovalsE_s;
        boost::mpi::gather(self->gcm_params.world, gcm_ovalsE_s,
            every_gcm_ovalsE_s, self->gcm_params.gcm_root);

        // Concatenate coupler inputs
        VectorMultivec gcm_ovalsE_s(concatenate(every_gcm_ovalsE_s));
#if 0
printf("BEGIN gcm_ovalsE_s nvar=%d\n", gcm_ovalsE_s.nvar);
for (size_t i=0; i<gcm_ovalsE_s.size(); ++i) {
    auto &iE_s(gcm_ovalsE_s.index[i]);
    double *vals(&gcm_ovalsE_s.vals[i*gcm_ovalsE_s.nvar]);
    printf("gcm_ovalsE_s[iE_s=%ld] =", iE_s);
    for (int j=0; j<gcm_ovalsE_s.nvar; ++j) printf(" %g", vals[j]);
    printf("\n");
}
printf("END gcm_ovalsE_s\n");
#endif

        // Couple on root!
        GCMInput out(
            self->couple(time_s,
                gcm_ovalsE_s, run_ice));

        // Split up the output (and 
        std::vector<GCMInput> every_outs(
            split_by_domain<DomainDecomposer_ModelE>(out,
                *self->domains, *self->domains));
#if 0
printf("BEGIN gcmce_couple_native() every_outs\n");
for (size_t i=0; i<every_outs.size(); ++i) {
    auto &gcmo(every_outs[i]);
    printf("    every_outs[%ld]: |gcm_ivalsA_s|=%ld, |gcm_ivalsE_s|=%ld\n", i,
        gcmo.gcm_ivalsAE_s[GridAE::A].size(),
        gcmo.gcm_ivalsAE_s[GridAE::E].size());
}
printf("END gcmce_couple_native() every_outs\n");
#endif

        // Scatter!
        boost::mpi::scatter(self->gcm_params.world, every_outs, out, self->gcm_params.gcm_root);


    } else {
        // =================== NOT MPI ROOT =============================
        // Send our input to root
        boost::mpi::gather(self->gcm_params.world, gcm_ovalsE_s, self->gcm_params.gcm_root);

        // Let root do the work...
        self->couple(time_s, gcm_ovalsE_s, run_ice);

        // Receive our output back from root
        boost::mpi::scatter(self->gcm_params.world, out, self->gcm_params.gcm_root);
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
    printf("BEGIN GCMCoupler_ModelE::update_gcm_ivals\n");
    auto nvar(out.nvar());    // A and E

    // Update gcm_ivalA variables...
    // Read from here...
    VectorMultivec const &gcm_ivalsA_s(out.gcm_ivalsAE_s[GridAE::A]);
    std::vector<std::unique_ptr<blitz::Array<double,2>>> &gcm_ivalsA(modele_inputs.gcm_ivalsA);
    if (gcm_ivalsA.size() != nvar[GridAE::A]) (*icebin_error)(-1,
        "gcm_ivalsA is wrong size: %ld vs. %d", gcm_ivalsA.size(), nvar[GridAE::A]);

    // Write to here...
    for (int iAE=0; iAE<GridAE::count; ++iAE) {
        if (nvar[iAE] < 0) (*icebin_error)(-1,
            "nvar[%d]=%d < 0, it should not be\n", iAE, nvar[iAE]);
    }

    for (int ivar=0; ivar<nvar[GridAE::A]; ++ivar) (*gcm_ivalsA[ivar]) = 0;
    for (size_t i=0; i<gcm_ivalsA_s.size(); ++i) {    // Iterate through elements of parallel arrays
        long iA = gcm_ivalsA_s.index[i];
        auto ij(gcm_regridder->indexing(GridAE::A).index_to_tuple<int,2>(iA));    // zero-based, alphabetical order
        int const i_f = ij[0]+1;    // C2F
        int const j_f = ij[1]+1;
        for (int ivar=0; ivar<nvar[GridAE::A]; ++ivar) {
            (*gcm_ivalsA[ivar])(i_f, j_f) +=
                gcm_ivalsA_s.vals[i*nvar[GridAE::A] + ivar];
        }
    }

    // Update gcm_ivalE variables...
    // Read from here...
    VectorMultivec const &gcm_ivalsE_s(out.gcm_ivalsAE_s[GridAE::E]);
    std::vector<std::unique_ptr<blitz::Array<double,3>>> &gcm_ivalsE(modele_inputs.gcm_ivalsE);
    if (gcm_ivalsE.size() != nvar[GridAE::E]) (*icebin_error)(-1,
        "gcm_ivalsE is wrong size: %ld vs. %d", gcm_ivalsE.size(), nvar[GridAE::E]);

    // Write to here...
    for (int ivar=0; ivar<nvar[GridAE::E]; ++ivar) (*gcm_ivalsE[ivar]) = 0;
    for (size_t i=0; i<gcm_ivalsE_s.size(); ++i) {
        long iE = gcm_ivalsE_s.index[i];
        auto ijk(gcm_regridder->indexing(GridAE::E).index_to_tuple<int,3>(iE));
        int const i_f = ijk[0]+1;    // C2F
        int const j_f = ijk[1]+1;
        int const ihc_ice = ijk[2];    // zero-based, just EC's known by ice model
        int const ihc_gcm_f = gcm_params.icebin_base_hc + ihc_ice + 1;    // 1-based, G

        for (int ivar=0; ivar<nvar[GridAE::E]; ++ivar) {
            (*gcm_ivalsE[ivar])(i_f, j_f, ihc_gcm_f) +=
                gcm_ivalsE_s.vals[i*nvar[GridAE::E] + ivar];
        }
    }
    printf("END GCMCoupler_ModelE::update_gcm_ivals\n");
}
// ============================================================================
// Update TOPO file during a coupled run



// NetCDF Logging Stuff





}}
