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
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <ibmisc/ncfile.hpp>
#include <ibmisc/string.hpp>
#include <ibmisc/f90blitz.hpp>
#include <ibmisc/math.hpp>
#include <icebin/modele/GCMCoupler_ModelE.hpp>
#include <icebin/modele/grids.hpp>
#include <icebin/contracts/contracts.hpp>
#include <icebin/domain_splitter.hpp>
#include <boost/filesystem.hpp>
#include <icebin/modele/GCMRegridder_ModelE.hpp>
#include <icebin/modele/hntr.hpp>
#include <spsparse/accum.hpp>
#include <spsparse/eigen.hpp>
#include <spsparse/SparseSet.hpp>

// See here to serialize objects with non-default constructor
//    http://www.boost.org/doc/libs/1_62_0/libs/serialization/doc/serialization.html#constructors
// http://www.boost.org/doc/libs/1_62_0/libs/serialization/doc/tutorial.html
// http://www.ocoudert.com/blog/2011/07/09/a-practical-guide-to-c-serialization/
// Binary archive that defines boost::archive::binary_oarchive
// and boost::archive::binary_iarchive
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>

using namespace std;
using namespace ibmisc;
using namespace netCDF;
using namespace spsparse;

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
void GCMCoupler_ModelE::_ncread(
    ibmisc::NcIO &ncio_config,
    std::string const &vname)        // comes from this->gcm_params
{
    GCMCoupler::_ncread(ncio_config, vname);

    auto config_info(get_or_add_var(ncio_config, vname + ".info", "int", {}));
    // Retrieve name of TOPO file (without Greenland, and on Ocean grid)
    get_or_put_att(config_info, ncio_config.rw, "topo_ocean", topoO_fname);

    // Replace the GCMRegridder with a wrapped version that understands
    // the ocean-vs-atmosphere grid complexity of ModelE
    GCMRegridder_ModelE *gcmA = new GCMRegridder_ModelE(gcm_regridder);
    gcm_regridder.reset(gcmA);

    // Allocate foceanOm0
//    Grid_LonLat *gridO = dynamic_cast<Grid_LonLat *>(&*gcmA->gcmO->gridA);
    GridSpec_LonLat *specO = dynamic_cast<GridSpec_LonLat *>(&*gcmA->gcmO->agridA.spec);
    foceanOm0.reference(blitz::Array<double,2>(specO->nlat(), specO->nlon()));
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
/** Reads from file ./config/icebin.nc */
extern "C"
GCMCoupler_ModelE *gcmce_new(
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
        auto config_info(get_or_add_var(ncio_config, "m.info", "int", {}));
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
        new blitz::Array<double,3>(f_to_c(var_f.to_blitz())));

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
        new blitz::Array<double,2>(f_to_c(var_f.to_blitz())));

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
        new blitz::Array<double,3>(f_to_c(var_f.to_blitz())));

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
    F90Array<double, 3> fhc,
    F90Array<int, 3> underice,
    F90Array<double, 3> elevE,
    F90Array<double, 2> focean,
    F90Array<double, 2> flake,
    F90Array<double, 2> fgrnd,
    F90Array<double, 2> fgice,
    F90Array<double, 2> zatmo)
{
    Topos *topos(&self->modele_inputs);
    topos->fhc.reference(f_to_c(fhc.to_blitz()));
    topos->underice.reference(f_to_c(underice.to_blitz()));
    topos->elevE.reference(f_to_c(elevE.to_blitz()));
    topos->focean.reference(f_to_c(focean.to_blitz()));
    topos->flake.reference(f_to_c(flake.to_blitz()));
    topos->fgrnd.reference(f_to_c(fgrnd.to_blitz()));
    topos->fgice.reference(f_to_c(fgice.to_blitz()));
    topos->zatmo.reference(f_to_c(zatmo.to_blitz()));
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
void gcmce_cold_start(GCMCoupler_ModelE *self, int yeari, int itimei, double dtsrc)
{
    printf("BEGIN gcmce_cold_start() yeari=%d, itimei=%d, dtsrc=%g\n", yeari, itimei, dtsrc);

    // This will cold-start initial conditions of the dynamic ice model
    self->dtsrc = dtsrc;

    // Call superclass cold_start()
    double const time_s = itimei * dtsrc;
    self->cold_start(
        ibmisc::Datetime(yeari,1,1), time_s);

    // b) Compute fhc, elevE
    // c) Compute ZATMO, FGICE, etc.
    self->update_topo(time_s, true);    // initial_timestep=true

    // d) Sync with dynamic ice model
    gcmce_couple_native(self, itimei, false);    // run_ice=false

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

    auto &indexingA(self->gcm_regridder->agridA.indexing);

    // domain uses alphabetical order, 0-based indexing...
    const auto base_hc(self->gcm_params.icebin_base_hc);
    const auto nhc_ice(self->gcm_regridder->nhc());

    auto &domainA(self->domainA);
printf("domainA size=%ld base_hc=%d  nhc_ice=%d\n", domainA.data.size(), base_hc, nhc_ice);

    for (int ihc=base_hc; ihc < base_hc + nhc_ice; ++ihc) {
        const int ihc_ice = ihc - base_hc;   // Use only IceBin HC's, zero-based indexing
        if (ihc_ice < 0) (*icebin_error)(-1,
            "ihc_ice cannot be <0: %d = %d - %d - 1\n", ihc_ice, ihc, base_hc);

        for (int j=domainA[1].begin; j < domainA[1].end; ++j) {
        for (int i=domainA[0].begin; i < domainA[0].end; ++i) {
            // i,j are 0-based indexes.
            if (self->modele_inputs.underice(ihc,j,i) != UI_ICEBIN) continue;
            long iE_s = self->gcm_regridder->indexingE.tuple_to_index(
                make_array(i, j, ihc_ice));

            if (iE_s < 0) (*ibmisc_error)(-1,
                "iE_s=%ld (from %d %d %d), it should not be negative\n", iE_s, i, j, ihc_ice);

            for (unsigned int ivar=0; ivar<self->gcm_outputsE.size(); ++ivar) {
                val[ivar] = (*self->modele_outputs.gcm_ovalsE[ivar])(ihc,j,i);
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
    self->update_topo(time_s, false);    // initial_timestep=false
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
    for (size_t ix=0; ix<gcm_ivalsA_s.size(); ++ix) {    // Iterate through elements of parallel arrays
        long iA = gcm_ivalsA_s.index[ix];
        auto ij(gcm_regridder->indexing(GridAE::A).index_to_tuple<int,2>(iA));    // zero-based, alphabetical order
        int const i = ij[0];
        int const j = ij[1];
        for (int ivar=0; ivar<nvar[GridAE::A]; ++ivar) {
            (*gcm_ivalsA[ivar])(j,i) +=
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
    for (size_t ix=0; ix<gcm_ivalsE_s.size(); ++ix) {
        long iE = gcm_ivalsE_s.index[ix];
        auto ijk(gcm_regridder->indexing(GridAE::E).index_to_tuple<int,3>(iE));
        int const i = ijk[0];
        int const j = ijk[1];
        int const ihc_ice = ijk[2];    // zero-based, just EC's known by ice model
        int const ihc_gcm = gcm_params.icebin_base_hc + ihc_ice;

        for (int ivar=0; ivar<nvar[GridAE::E]; ++ivar) {
            (*gcm_ivalsE[ivar])(ihc_gcm,j,i) +=
                gcm_ivalsE_s.vals[i*nvar[GridAE::E] + ivar];
        }
    }
    printf("END GCMCoupler_ModelE::update_gcm_ivals\n");
}
// ============================================================================
// Update TOPO file during a coupled run



/** Adds ice sheet information to an FOCEAN read from Gary's topo files.

@param foceanOp foceanO as read in from TOPO files, with ice sheets
    removed.  Starting out, foceanOp should be 0 or 1 everywhere.
    This function will change foceanOp in areas of ice sheets, setting
    to values in the range [0,1]
*/
TODO: Remove update_foceanOp
static void update_foceanOp(
GCMRegridder *gcmO,
std::vector<ElevMask<1>> const &elevmasks,
blitz::Array<double,1> &foceanOp,    // OUT: 0-based array
blitz::Array<char,1> &changedO)    // OUT
{
//    GCMRegridder *gcmA = &*gcmc->gcm_regridder;

This IS MUCH MORE EASILY DONE as a byproduct of merged EOpvAOp (without global ice), which must be computed anyway.

    auto nO(gcmO->nA());

    // --------------------- Compute fcontOp_d (and foceanOp)
    for (size_t sheet_index=0; sheet_index < gcmO->ice_regridders().size(); ++sheet_index) {

        // Construct an elevmaskI for CONTINENTAL land (not just ice sheet)
        // (==elevI on continent, nan on ocean)
        auto &emI(elevmasks[sheet_index]);
        int nI = emI.elev.extent(0);
        blitz::Array<double,1> elevmaskI(nI);
        for (int iI=0; iI<nI; ++iI) {
            auto const m(emI.mask(iI));
            elevmaskI(iI) = (m == IceMask::ICE_FREE_OCEAN ? nan : emI.elev(iI));
        }

        // Get OvI for continental cells
        std::unique_ptr<RegridMatrices_Dynamic> rmO(gcmO->regrid_matrices(sheet_index, elevmaskI));
        SparseSetT dimO, dimI;
        RegridParams paramsO(paramsA);
            paramsO.scale = false;
            paramsO.correctA = false;
        auto OvI(rmO->matrix_d("AvI", {&dimO, &dimI}, paramsO));

        // Don't need to set up the mask on I ourselves; this is already
        // built into the OvI matrix.  The mask, taken from PISM, includes
        // all bare land and ice-covered areas.
        // See: pygiss/giss/pism.py   _get_landmask()
        //    (used by write_icebin_in_base.py)
        blitz::Array<double,1> fcontI_d(dimI.dense_extent());
        fcontI_d = 1.0;

        // Compute fcontOp (for this ice sheet only)
        TmpAlloc tmp;

        blitz::Array<double,1> fcontOp_d(OvI->apply(fcontI_d, 0., true, tmp));    // force_conservation set to true by default, and it probably doesn't matter; but I think it should be false here.
TODO: This will produce AREA of continent, not FRACTION of continent

        // Interpolate into foceanOp_s
        IceRegridder *iceO = &*gcmO->ice_regridders()[sheet_index];
        ibmisc::Proj_LL2XY proj(iceO->agridI.sproj);

        GridSpec_LonLat const *gridO = dynamic_cast<GridSpec_LonLat *>(&*gcmO->agridA.spec);
        for (int iO_d=0; iO_d<fcontOp_d.extent(0); ++iO_d) {



            // fcont will be 0 or 1; must multiply by fraction of gridcell covered by continent.
            auto const iO_s = dimO.to_sparse(iO_d);
            // double const area = gridO->cells.at(iO_s)->proj_area(&proj);
            double const area = iceO->gridA_proj_area(iO_d);

            foceanOp(iO_s) = foceanOp(iO_s) - round_mantissa(fcontOp_d(iO_d) / area, 12);
            changedO(iO_s) = 1;    // true
        }

#if 0
// Write debugging NetCDF file
{NcIO ncio("x.nc", 'w');

    Grid_LonLat const *gridO(cast_Grid_LonLat(&*gcmO->gridA));
    auto shapeO(blitz::shape(gridO->nlat(), gridO->nlon()));


    // fcontOp
    blitz::Array<double,1> fcontOp(nO);
    fcontOp = nan;
    for (int i=0; i<dimO.dense_extent(); ++i) {
        fcontOp(dimO.to_sparse(i)) = fcontOp_d(i);
    }

    // (IM, IM)
    auto dimsO(get_or_add_dims(ncio, {"jm", "im"}, {gridO->nlat(), gridO->nlon()}));
    auto fcontOp2(reshape<double,1,2>(fcontOp, shapeO));
    auto foceanOp2(reshape<double,1,2>(foceanOp, shapeO));

    ncio_blitz(ncio, fcontOp2, false, "fcontOpx", dimsO);
    ncio_blitz(ncio, foceanOp2, false, "foceanOp", dimsO);

    ncio.close();
}
#endif

    }


}

/** Adds ice sheets to Gary's TOPO file */
static void update_fgiceO(
GCMRegridder *gcmO,
std::vector<ElevMask<1>> const &elevmasks,
blitz::Array<double,1> &fgiceO,    // OUT: 0-based array
blitz::Array<char,1> &changedO)    // OUT
{
    auto nO(gcmO->nA());

    // --------------------- Compute fgiceO
    for (size_t sheet_index=0; sheet_index < gcmO->ice_regridders().size(); ++sheet_index) {
        TmpAlloc tmp;
    
        // Construct an elevmaskI for CONTINENTAL land and ICE SHEET
        // (==elevI on continent, nan on ocean)
        auto &emI(elevmasks[sheet_index]);
        int nI = emI.elev.extent(0);
        blitz::Array<double,1> elevmaskI(nI);
        for (int iI=0; iI<nI; ++iI) {
            auto const m(emI.mask(iI));
            elevmaskI(iI) = (
                m==IceMask::GROUNDED_ICE || m==IceMask::FLOATING_ICE || m==IceMask::BARE_GROUND ?
                emI.elev(iI) : nan);
        }

        // Get OvI for ice cells
        std::unique_ptr<RegridMatrices_Dynamic> rmO(gcmO->regrid_matrices(sheet_index, elevmaskI));
        SparseSetT dimO, dimI;
        RegridParams paramsO;
            paramsO.scale = true;
            paramsO.correctA = false;
        auto OvI(rmO->matrix_d("AvI", {&dimO, &dimI}, paramsO));

        // Don't need to set up the mask on I ourselves; this is already
        // built into the OvI matrix.  The mask, taken from PISM, includes
        // all bare land and ice-covered areas.
        // See: pygiss/giss/pism.py   _get_landmask()
        //    (used by write_icebin_in_base.py)
        blitz::Array<double,1> fgiceI_d(dimI.dense_extent());
        fgiceI_d = 1.0;

        // Compute fgiceO (for this ice sheet only)
        blitz::Array<double,1> fgiceO_d(OvI->apply(fgiceI_d, 0., true, tmp));    // force_conservation set to true by default, and it probably doesn't matter; but I think it should be false here.

        // Interpolate into global FGICE
        for (int iO_d=0; iO_d<fgiceO_d.extent(0); ++iO_d) {
            auto const iO_s = dimO.to_sparse(iO_d);
            fgiceO(iO_s) += fgiceO_d(iO_d);        // Will have some rounding error on #s that should ==1.0
            changedO(iO_s) = 1;    // true
        }

    }
}

/** Construct merged data on Ocean grid; ready as input for make_tooa
(variables named after make_topoa.cpp) */
class MergeTopoO {
    EigenSparseMatrixT AOvEO;
    HnterSpec hspec;
    blitz::Array<double,2> fgiceOp;
    blitz::Array<double,2> fgiceOm;
    blitz::Array<double,3> elevEO;
};

static std::vector<std::string> const otopo_vars
    {"focean", "flake", "fgrnd", "fgice", "zatmo", "hlake", "zicetop"};
std::map<std::string, std::unique_ptr<blitz::Array<double,2>>> varsO;


get_OvI(
GCMRegridder *gcmO,
RegridParams const &paramsA,
int sheet_index,
ElevMask<1> &emI,
SparseSetT &dimOp,
bool include_ice,
bool include_bedrock)
{
        // Construct an elevmaskI for ice only
        // (==elevI on ice, nan on ocean or bare land)
        auto &emI(elevmasks[sheet_index]);
        int nI = emI.elev.extent(0);
        blitz::Array<double,1> elevmaskI(nI);
        for (int iI=0; iI<nI; ++iI) {
            switch(emI.mask(iI)) {
                case IceMask::GROUNDED_ICE :
                case IceMask::FLOATING_ICE :
                    elevmaskI(iI) = (include_ice ? emI.elev(iI) : nan);
                break;
                case IceMask::ICE_FREE_OCEAN :
                case IceMask::UNKNOWN :
                    elevmaskI(iI) = nan;
                break;
                case IceMask::ICE_FREE_BEDROCK :
                    elevmask(iI) = (include_bedrock ? emI.elev(iI) : nan);
                break;
            }

        }


        SparseSet dimOp,dimI;
        std::unique_ptr<RegridMatrices_Dynamic> rmO(
            gcmO->regrid_matrices(sheet_index, elevmaskI));
        RegridParams paramsO(paramsA);
            paramsO.scale = false;
            paramsO.correctA = false;
       return rmO->matrix_d("AvI", {&dimOp, &dimI}, paramsO);
}



void merge_topoO(
// TOPOO arrays to merge into (global ice)
blitz::Array<double,2> const &areaO,    // Area of each ocean gridcell
blitz::Array<double,2> &foceanOp,    // Fractional FOCEAN
blitz::Array<double,2> &foceanOm,     // Rounded FOCEAN
blitz::Array<double,2> &flakeO,
blitz::Array<double,2> &fgrndO,
blitz::Array<double,2> &fgiceO,

blitz::Array<double,2> &zatmoO,
blitz::Array<double,2> &hlakeO,
blitz::Array<double,2> &zicetopO,
// Local ice to merge in...
GCMRegridder *gcmO,
std::vector<ElevMask<1>> const &elevmasks,
double const eq_rad,    // Radius of the earth
std::vector<ElevMask<1>> const elevmasks,    // elevation and cover types for each ice sheet


//EigenSparseMatrixT const &AOmvEOm,
//SparseSetT const &dimAOm,
//SparseSetT const &dimEOm)
{
    SparseSetT dimEOp, dimAOp;
    auto nO = gcmO->nA();

    // -------------- Multiply key variables by (1-foceanOp), to accumulate later
    auto &zicetopO_area(zicetopO);
    for (long iO_s=0; iO_s < nO; ++iO) {
        zatmoO_frac(iO_s) *= (1. - foceanOp(iO_s));
        zicetopO_frac(iO_s) *= fgiceOp(iO_s);
    }

#if 0
    // ------------------ Update land surface fractions from local ice sheets
    fcontO_local = 0.;
    {
        // Get area of non-ocean (fcont/fgrnd) in each gridcell we touch
        SparseSetT dimEOp, dimAOp;
        EigenMatrixT EOpvAOp(compute_EOpvAOp_merged(
            {&dimEOp, &dimAOp},
            paramsA, gcmO, eq_rad, elevmasks,
            false, true, true));  // local ice only, include bare ground
        auto area_contO(sum(EOpvAOp_unscaled, 1, '+'));

        // Update foceanOp
        for (size_t iAO_d=0; iAO_d < dimAOp.size(); ++iAO_d) {
            long const iAO_s = dimAO.to_sparse(iAO_d);

            fcontOp_local(iAO_d) += area_contO(iAO_d) / areaO(iAO_s);
//            double fcontOp_old = 1.0 - foceanOp(iAO_s);
//            double fcontOp_local = area_contO(iAO_d) / areaO(iAO_s);
//            foceanOp(iAO_s) = 1.0 - (fcontOp_old + fcontOp_local);


        }
    }
#endif


    // -------------- Compute ZATMO (elevO), etc.
    for (size_t sheet_index=0; sheet_index < gcmO->ice_regridders().index.size(); ++sheet_index) {

        // Construct elevI in dense indexing space
        ElevMask<1> const &emI(elevmasks[sheet_index]);
        blitz::Array<double,1> elevI_d(dimI.size());
        for (size_t iI_d=0; iI_d < dimIp.size(); ++iI_d) {
            auto iI_s = dimI.to_sparse(iI_d);
            elevI_d(iO_d) = emI.elev(iI_s);
        }

        blitz::Array<double,1> dacontO(nO);    // increase in area of continent [m^2]
        dacontO = 0.;
        blitz::Array<double,1> dagiceO(nO);     // increase in area of landice [m^2]
        dagiceO = 0.;

        // Update zicetopO from ice-only coverage
        {TmpAlloc tmp;
            SparseSetT dimOp;
            auto OvI_ice(get_OvI(gcmO, paramsA, sheet_index, emI, dimOp, true, false); // unscaled

            for (size_t iO_d=0; iO_d < dimOp.size(); ++iO_d) {
                auto iO_s = dimO.to_sparse(iO_d);
            auto elevO_d(OvI_ice->apply(elevI_d, nan, true, tmp));
            dagiceO(iO_s) += OvI->wM(iO_d);
            zicetopO_frac(iO_s) += elevO_d(iO_d) * OvI->wM(iO_d) / area(iO_s);;
        }

        {TmpAlloc tmp;
            SparseSetT dimOp;
            auto OvI_ice(get_OvI(gcmO, paramsA, sheet_index, emI, dimOp, true, false); // unscaled

            for (size_t iO_d=0; iO_d < dimOp.size(); ++iO_d) {
                auto iO_s = dimO.to_sparse(iO_d);
            auto elevO_d(OvI_ice->apply(elevI_d, nan, true, tmp));
            dacontO(iO_s) += OvI->wM(iO_d);
            zatmoO_frac(iO_s) += elevO_d(iO_d) * OvI->wM(iO_d) / area(iO_s);;
        }





            SparseSetT dimOp_cont;
            auto OvI_cont(get_OvI(gcmO, paramsA, sheet_index, emI, dimOp_cont, true, true);

            // Compute elevO (dense indexing)
            auto elevO_cont_d(OvI_cont->apply(elevI_d, nan, true, tmp));

        for (size_t iO_d=0; iO_d < dimOp.size(); ++iO_d) {
            auto iO_s = dimO.to_sparse(iO_d);
            // Replace ocean in original (global) ZATMO with
            // (ice+land) from local ice sheet.
            zicetopO_area(iO_s) += elevO_d(iO_d) * OvI->wM(iO_d);
        }
    }


    // Create a rounding plan
    int8_t const TO_NONE = 0;
    int8_t const TO_OCEAN = 1;
    int8_t const TO_LAND = 2;
    blitz::Array<int8_t,2> SET_TO(IM,JM);
    SET_TO = TO_NONE;


    for (size_t iAO_d=0; iAO_d < dimAOp.size(); ++iAO_d) {
        long const iAO_s = dimAO.to_sparse(iAO_d);

        // Adjust fgrnd to make it all sum to 1; and round foceanOm at the same time
        flakeO(iAO_s) = 0.;
        double const focean = foceanOp(iAO_s);
        if (focean > 0. && focean < 0.5) {
            SET_TO(iAO_s) = TO_LAND;
        } else if (focean >= 0.5 && focean < 1.0) {
            SET_TO(iAO_s) = TO_OCEAN;
        }
    }


    // ---------------- Eliminate single-cell oceans
    GridSpec_LonLat const &specO(cast_GridSpec_LonLat(*gcmO->agridA.spec));
    auto shapeO(blitz::shape(specO.nlat(), specO.nlon()));
    auto foceanOm2(reshape<double,1,2>(foceanOm, shapeO));

    auto const im(specO.nlon());
    auto const jm(specO.nlat());
    std::array<long,2> ijO;
    for (int j=0; j<JM; ++j) {
    for (int i=0; i<IM; ++i) {

        int iO(gcmO->agridA.indexing.tuple_to_index(ijO));
        auto const i(ijO[0]);
        auto const j(ijO[1]);

        // Avoid edges, where indexing is more complex (and we don't need to correct anyway)
        if ((i==0) || (i==im-1) || (j==0) || (j==jm-1)) continue;

        if (((foceanOp2(j-1,i) == 0.) || (SET_TO2(j-1,i)==TO_LAND))
            && ((foceanOp2(j+1,i) == 0.) || (SET_TO2(j+1,i)==TO_LAND))
            && ((foceanOp2(j,i-1) == 0.) || (SET_TO2(j,i-1)==TO_LAND))
            && ((foceanOp2(j,i+1) == 0.) || (SET_TO2(j,i+1)==TO_LAND)))
        {
            // This cell is surrounded by land
            if (foceanOp2(j,i) == 0.) SET_TO2(j,i)=TO_NONE;   // defensive: make sure we don't change it
            else SET_TO2(j,i) = TO_LAND;
        }
    }

    // Apply the rounding plan
    for (int iAO_s=0; iAO_s < nAO; ++iAO_s) {
        // Adjust fgrnd to make it all sum to 1; and round foceanOm at the same time
        switch(SET_TO(iAO_s)) {
            case TO_OCEAN :
                foceanOm(iAO_s) = 1.;
                fgiceO(iAO_s) = 0.;
                fgrndO(iAO_s) = 0.;
                zatmoO(iAO_s) = 0.;    // All ocean!
                zicetopO(iAO_s) = 0.;    // Doesn't matter, there's no ice
                hlakeO(iAO_s) = 0.;    // All ocean at sea level
            break;
            case TO_LAND :
                foceanOm(iAO_s) = 0.;
                fgiceO(iAO_s) = round_mantissa(fgiceO(iAO_s), 3);    // will be ==1.0 when needed
                fgrndO(iAO_s) = 1. - fgiceO(iAO_s);    // Should have pure zeros, not +-1e-17 stuff
                // zicetopO doesn't change
                // zatmoO: Adjust for land now filling entire cell
                zatmoO(iAO_s) /= (1. - foceanOm(iAO_s));
            break;
        }
    }

}

void make_topoA(
GCMRegridder_ModelE *gcmA,    // Gets updated with new fcoeanOp, foceanOm
blitz::Array<double,2> const &foceanOm,     // Rounded FOCEAN
blitz::Array<double,2> const &flakeO,
blitz::Array<double,2> const &fgrndO,
blitz::Array<double,2> const &fgiceO,
blitz::Array<double,2> const &zatmoO,
//SparseSetT const &dimO,    // Tells us which grid cells in O were changed.
blitz::Array<double,2> const &foceanA,    // Rounded FOCEAN
blitz::Array<double,2> const &flakeA,
blitz::Array<double,2> const &fgrndA,
blitz::Array<double,2> const &fgiceA,
blitz::Array<double,2> const &zatmoA,

)
{

    // =====================================================
    // Regrid TOPO to Atmosphere grid
    HntrSpec const &hntrA(cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr);
    HntrSpec const &hntrO(cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr);
    Hntr hntrAvO(17.17, hntrA, hntrO);

    TupleListT<2> AvO_tp;
    hntrAvO.scaled_regrid_matrix(spsparse::accum::ref(AvO_tp));
    EigenSparseMatrixT AvO_e(hntrA.size(), hntrO.size());
    AvO_e.setFromTriplets(AvO_tp.begin(), AvO_tp.end());

    fgiceA = 0;

    map_eigen_colvector(foceanA) = AvO_e * map_eigen_colvector(foceanOm);
    map_eigen_colvector(flakeA) = AvO_e * map_eigen_colvector(flakeO);
    map_eigen_colvector(fgrndA) = AvO_e * map_eigen_colvector(fgrndO);
    map_eigen_colvector(fgiceA) = AvO_e * map_eigen_colvector(fgiceO);
    map_eigen_colvector(zatmoA) = AvO_e * map_eigen_colvector(zatmoO);
    map_eigen_colvector(elevA) = AvO_e * map_eigen_colvector(elevO);


    TupleListT<2> AvE_global_tp;
    blitz::Array<double,1> elevE_global(nE);


}
// ======================================================================

THE PROBLEM:
   Re-running make_topoa requires re-building AvE mismatched matrix for 1-minute ice cover, which is slow and memory intensive.  Ways around this problem:

1. make_topoa only needs mismatched AvE, which can be generated from non-mismatched matrices on ocean grid.  Inspection of GCMRegridder_ModelE indicates that computing AvE on atmosphere grid (compute_AAmvEAm) requires EvA from ocean grid (EOpvAOp) in EOmvAOm_Unscaled().  It also requires wAOm (PISM view of ocean grid), whose calculation nominally requires AvI on ocean grid.  HOWEVER... this is just used to compute wAOp (ice model view on ocean grid), which can be pre-computed and stored, without further requiring large ice-size vectors.

wAOp from global ice can be merged with wAOp from ice sheets (add them together), resulting in a global wAOp including ice sheets.  Which can then be used to make TOPOA.  (NOTE: wAOp might be really similar to FOCEANp from TOPOO file; in which case, maybe it is already stored.)

So... the general approach is to:

A. Load pre-computed global ice stuff on ocean grid:
 1. Load TOPOO
 2. load global_ecO.nc files (global EC matrices on ocean grid), which contains OAvOE.
 3. Load wOp (wAOp) from global_ecA_mm.nc (of which it is a byproduce).  Or decide it can be recomputed from foceanOp (generated in make_topoo.cpp)

B. Internally generate TOPO and regrid matrix stuff for the ice sheets.

C. Merge (B) into (A)

D. Re-run make_topoa to finish the job.  It will be supplied our own (merged) wOp and AvE, rather than computing it from global_ecA_mm.



Overall merge procedure will require:
  * topoo.nc
  * topoa.nc

Should we re-run (essentially) make_topoa with a merge involved?  Or should we do something similar but different in code here?  Probably the latter, since we are ADDING to the result of make_topoa.nc
ANSWER: We should re-run (code factored out of) make_topoa.cpp.  This will ensure:
  a) Changes to TOPOO propagate through correctly.
  b) A single codebase is used for dynamic (internal) vs. static production of TOPO files.
For this, make_topoa will need to be augmented to accept per-ice sheet inputs.
  (i.e. accept a merged AvE for all the ice sheets, and merge that into the global AvE).

HOWEVER... global_ec_mm requires etopo1 and TOPOO, so it is a real bear.  And it is needed for make_topoa.cpp.  So 

void update_topo(
    // ====== INPUT parameters
    GCMRegridder_ModelE *gcmA,    // Gets updated with new fcoeanOp, foceanOm
    std::string const &topoO_fname,    // Name of Ocean-based TOPO file (aka Gary)
    std::vector<ElevMask<1>> const &elevmasks,
    std::vector<std::array<double,3>> const &sigmas,
    bool initial_timestep,    // true if this is the first (initialization) timestep
    std::vector<HCSegmentData> const &hc_segments,
    std::string const &primary_segment,    // Normally "ec"
    // ===== OUTPUT parameters (variables come from GCMCoupler); must be pre-allocated
    Topos &topoA,
    blitz::Array<double,2> foceanOm0)
{    // BEGIN update_topo

    if (!initial_timestep) (*icebin_error)(-1,
        "GCMCoupler_ModelE::update_topo() currently only works for the initial call");

    HCSegmentData const &legacy(get_segment(hc_segments, "legacy"));
    HCSegmentData const &sealand(get_segment(hc_segments, "sealand"));
    HCSegmentData const &ec(get_segment(hc_segments, "ec"));

    GCMRegridder *gcmO = &*gcmA->gcmO;
    auto nA = gcmA->nA();
    auto nE = gcmA->nE();
    auto nO = gcmO->nA();
    auto nhc_ice = gcmA->nhc();
    int nhc_gcm = ec.base + nhc_ice;

    // Convert TOPO arrays to 1-D zero-based indexing
    // ...on elevation grid
    blitz::TinyVector<int,2> const shape_E2(nhc_gcm, nA);
    blitz::Array<double,2> fhcE2(reshape(topoA.fhc, shape_E2));
    blitz::Array<int,2> undericeE2(reshape(topoA.underice, shape_E2));
    blitz::Array<double,2>  elevE2(reshape(topoA.elevE, shape_E2));
    // ...on atmosphere grid
    auto foceanA(reshape1(topoA.focean));
    auto flakeA(reshape1(topoA.flake));
    auto fgrndA(reshape1(topoA.fgrnd));
    auto fgiceA(reshape1(topoA.fgice));
    auto zatmoA(reshape1(topoA.zatmo));

    // Read the original topo file [Ocean grid]
    NcIO ncio(topoO_fname, 'r');
    auto foceanO2(nc_read_blitz<double,2>(ncio.nc, "FOCEAN"));
    auto flakeO2(nc_read_blitz<double,2>(ncio.nc, "FLAKE"));
    auto fgrndO2(nc_read_blitz<double,2>(ncio.nc, "FGRND"));
    auto fgiceO2(nc_read_blitz<double,2>(ncio.nc, "FGICE"));
    auto zatmoO2(nc_read_blitz<double,2>(ncio.nc, "ZATMO"));

    blitz::Array<double,1> foceanO(reshape1(foceanO2));
    blitz::Array<double,1> flakeO(reshape1(flakeO2));
    blitz::Array<double,1> fgrndO(reshape1(fgrndO2));
    blitz::Array<double,1> fgiceO(reshape1(fgiceO2));
    blitz::Array<double,1> zatmoO(reshape1(zatmoO2));
    ncio.close();

    // Keep track of which gridcells have been changed
    blitz::Array<char,1> changedO(nO);
    changedO = 0;

    // --------------------------------------
    // Add ice sheet to foceanO (and call it foceanOp)
    blitz::Array<double,1> &foceanOp(gcmA->foceanAOp);
    blitz::Array<double,1> &foceanOm(gcmA->foceanAOm);
    foceanOp = foceanO;
    update_foceanOp(gcmO, elevmasks, foceanOp, changedO);

    // --------------------------------------
    // Add ice to the surface type
    update_fgiceO(gcmO, elevmasks, fgiceO, changedO);

    // --------------------------------------
    // Adjust fgrnd to make it all sum to 1; and round foceanOm at the same time
    for (int i=0; i<nO; ++i) {
        if (changedO(i)) {
            flakeO(i) = 0.;
            if (foceanOp(i) >= 0.5) {
                foceanOm(i) = 1.;
                fgiceO(i) = 0.;
                fgrndO(i) = 0.;
            } else {
                foceanOm(i) = 0.;
                fgiceO(i) = round_mantissa(fgiceO(i), 3);    // will be ==1.0 when needed
                fgrndO(i) = 1. - fgiceO(i);    // Should have pure zeros, not +-1e-17 stuff
            }
        } else {
            foceanOm(i) = foceanOp(i);
        }
    }

    // ----------------------------------------------
    // Eliminate single-cell oceans
    GridSpec_LonLat const &specO(cast_GridSpec_LonLat(*gcmO->agridA.spec));
    auto shapeO(blitz::shape(specO.nlat(), specO.nlon()));
    auto foceanOm2(reshape<double,1,2>(foceanOm, shapeO));

    auto const im(specO.nlon());
    auto const jm(specO.nlat());
    // Avoid edges, where indexing is more complex (and we don't need to correct anyway)
    std::array<int,2> ijO{1,1};
    for (ijO[1]=1; ijO[1]<jm-1; ++ijO[1]) {
    for (ijO[0]=1; ijO[0]<im-1; ++ijO[0]) {
        int iO(gcmO->agridA.indexing.tuple_to_index(ijO));
        auto const i(ijO[0]);
        auto const j(ijO[1]);

        if (changedO(iO)) {
            if (foceanOm2(j,i) == 1. && foceanOm2(j-1,i)==0. && foceanOm2(j,i-1) == 0. && foceanOm2(j+1,i) == 0. & foceanOm2(j,i+1) == 0) {
                // Repeat of if-body above
                foceanOm(iO) = 0.;
                fgiceO(iO) = round_mantissa(fgiceO(iO), 3);
                fgrndO(iO) = 1. - fgiceO(iO);
            }
        }
    }}


    // Store the initial FOCEAN for ModelE, since it cannot change later.
    if (initial_timestep) {
        auto foceanOm0_1(reshape1(foceanOm0));
        foceanOm0_1 = foceanOm;
    }



# if 0
// Debugging NetCDF file
{NcIO ncio("y.nc", 'w');

    GridSpec_LonLat const &specA(cast_GridSpec_LonLat(*gcmA->agridA.spec));
    auto shapeA(blitz::shape(specA.nlat(), specA.nlon()));
    GridSpec_LonLat const &specO(cast_GridSpec_LonLat(*gcmO->agridA.spec));
    auto shapeO(blitz::shape(specO.nlat(), specO.nlon()));

    // (IM, IM)
    auto dimsO(get_or_add_dims(ncio, {"jmO", "imO"}, {specO.nlat(), specO.nlon()}));
    auto dimsA(get_or_add_dims(ncio, {"jmA", "imA"}, {specA.nlat(), specA.nlon()}));

    auto foceanOp2(reshape<double,1,2>(foceanOp, shapeO));
    ncio_blitz(ncio, foceanOp2, false, "foceanOp", dimsO);

    auto foceanOm2(reshape<double,1,2>(foceanOm, shapeO));
    ncio_blitz(ncio, foceanOm2, false, "foceanOm", dimsO);

    ncio.close();
}
#endif


====================== Up to here should be pretty straightforward.
Put it into a merge_topoo, whose output can be directly comparable to make_topoo
Output is a merged TOPOO (no FHC)



    // ----------------------------------------------------------
    // ----------------------------------------------------------
    // Now we are ready to use regrid matrices

    // =====================================================
    // Regrid TOPO to Atmosphere grid
    HntrSpec const &hntrA(cast_GridSpec_LonLat(*gcmA->agridA.spec).hntr);
    HntrSpec const &hntrO(cast_GridSpec_LonLat(*gcmA->gcmO->agridA.spec).hntr);
    Hntr hntrAvO(17.17, hntrA, hntrO);

    TupleListT<2> AvO_tp;
    hntrAvO.scaled_regrid_matrix(spsparse::accum::ref(AvO_tp));
    EigenSparseMatrixT AvO_e(hntrA.size(), hntrO.size());
    AvO_e.setFromTriplets(AvO_tp.begin(), AvO_tp.end());

    fgiceA = 0;

    map_eigen_colvector(foceanA) = AvO_e * map_eigen_colvector(foceanOm);
    map_eigen_colvector(flakeA) = AvO_e * map_eigen_colvector(flakeO);
    map_eigen_colvector(fgrndA) = AvO_e * map_eigen_colvector(fgrndO);
    map_eigen_colvector(fgiceA) = AvO_e * map_eigen_colvector(fgiceO);
    map_eigen_colvector(zatmoA) = AvO_e * map_eigen_colvector(zatmoO);

    TupleListT<2> AvE_global_tp;
    blitz::Array<double,1> elevE_global(nE);


Before merging in AvE from ice sheets, we need to start with AvE from topoa.nc

------------------------------
This code is just for ice sheets, not global ice (which has super-big EvI that never changes; so we might as well use the pre-computed fhc):
    // Compute elevE and AvE (aka fhc)
    SparseSetT dimA_global;
    for (size_t sheet_index=0; sheet_index < gcmO->ice_regridders().index.size(); ++sheet_index) {
        TmpAlloc tmp;

        ElevMask<1> &emI(elevmasks[sheet_index]);
        int nI(emI.elev.extent(0));

        // Construct an elevmaskI for ice sheet, =nan off ice sheet
        blitz::Array<double,1> elevmaskI(nI);
        for (int iI=0; iI<nI; ++iI) {
            auto const m(emI.mask(iI));
            if (m==IceMask::GROUNDED_ICE || m==IceMask::FLOATING_ICE) {
                elevmaskI(iI) = emI.elev(iI);
            } else {
                elevmaskI(iI) = nan;
            }
        }

        // Get regrid matrice needed to compute global stuff
        RegridParams params;
            params.scale = true;
            params.correctA = false;
            params.sigma = sigmas[sheet_index];    // TODO: Set smoothing!
        std::unique_ptr<RegridMatrices_Dynamic> rmA(gcmA->regrid_matrices(sheet_index, elevmaskI, params));
        SparseSetT dimA, dimE, dimI;
//        auto AvI(rmA->matrix_d("AvI", {&dimA, &dimI}, params));   NOT USED
//        auto EvI(rmA->matrix_d("EvI", {&dimE, &dimI}, params));   NOT USED
        auto AvE(rmA->matrix_d("AvE", {&dimA, &dimE}, params));

        // Merge local (per-ice sheet) into global AvE
TODO (BUG FIX): Merging like this requires UNSCALED matrices!
        spsparse::spcopy(
            spsparse::accum::to_sparse(std::array<SparseSetT *,2>{&dimA, &dimE},
            spsparse::accum::ref(AvE_global_tp)),
            *AvE->M, true);


        // Densify elevmaskI
        blitz::Array<double,1> elevmaskI_d(dimI.dense_extent());
        for (int i_d=0; i_d<dimI.dense_extent(); ++i_d) {
            elevmaskI_d(i_d) = elevmaskI(dimI.to_sparse(i_d));
        }

        // elevE is the mean elevation of each elevation class's actualy ice.
        // It will be used to create elevA (global), which is used for ZATMO, etc.
        auto elevE_d(EvI->apply(elevmaskI_d, nan, true, tmp));

There is no point in this elevE variable.  elevE -> elevE_global -> elevA_global.  So we might as well just compute elevA_global directly, using AvI instead of EvI.  We need elevA for ZATMO, etc.  THe elevE array should be set to the hcdefs points, which are known beforehand.  AvI is smaller than EvI.


        // Merge local and global elevE
        spsparse::spcopy(
            spsparse::accum::to_sparse(std::array<SparseSetT *,1>{&dimE},
            spsparse::accum::blitz_existing(elevE_global)),
            elevE_d, true);

        // Add to dimA_global
        for (int iA_d=0; iA_d<dimA.dense_extent(); ++iA_d) {
            int iA_s = dimA.to_sparse(iA_d);
            dimA_global.add_dense(iA_s);
        }
    }
------------------------------


    // Create matrix that works directly on sparse-indexed vectors
    EigenSparseMatrixT AvE_global_e(nA, nE);
    AvE_global_e.setFromTriplets(AvE_global_tp.begin(), AvE_global_tp.end());
    EigenColVectorT elevA_global_e(AvE_global_e * map_eigen_colvector(elevE_global));
    auto elevA_global(to_blitz(elevA_global_e));    // Sparse indexing

    // =======================================================
    // ----------- Set TOPO variables

    // =======================================================
    // ----------- Set up elevation class structure

Here we need to rearrange how legacy is computed, based on make_topoa.cpp.
And we need to add separate EC segments for global ECs and per-ice-sheet ECs.
(Yes we can compress in the future if we decide we have too many ECs)

So segments are: legacy,sealand,global,icesheets

In summary, which matrices does ModelE+PISM need?
EvA: trivial
AvE: ---> fhc
EvI,IvE: Only for ice sheets
AvI,IvA: Do we need this?  If so, it is only for ice sheets

This code should be analogous to make_topoa.cpp; exepct it adds ice sheets to what it otherwise computers (fhc,elevE,undericeE).

    // Set up elevation class segments: fhc, underice, elevmaskI
    fhcE2 = 0;
    undericeE2 = 0;
    elevE2 = 0;

    double const zero = 1.e-30;    // Non-zero, yet adds nothing
    double const legacy_mult = (primary_segment == "legacy" ? 1.0 : zero);
    double const sealand_mult = (primary_segment == "sealand" ? 1.0 : zero);
    double const ec_mult = (primary_segment == "ec" ? 1.0 : zero);

    // ------- Segment 0: Legacy Segment
    for (int ihc=legacy.base; ihc<legacy.base+legacy.size; ++ihc) {
        // Full domain
        for (int iA_s=0; iA_s<nA; ++iA_s) {
            if (fgiceA(iA_s) != 0) {
                fhcE2(ihc,iA_s) = 1.0;
                undericeE2(ihc,iA_s) = UI_NOTHING;
                elevE2(ihc,iA_s) = zatmoA(iA_s);
            }
        }

        // overlay...
        for (int iA_d=0; iA_d<dimA_global.dense_extent(); ++iA_d) {
            int iA_s = dimA_global.to_sparse(iA_d);

This is wrong for elevE2.  It must be elevA_global(iA_s) * fhc
            elevE2(ihc,iA_s) = elevA_global(iA_s);
            fhcE2(ihc,iA_s) = legacy_mult;    // Legacy ice for Greenland and Antarctica.
        }



    }

    // ------- Segment 1: SeaLand Segment
    // ihc=0: Non-ice portion of grid cell at sea level

    // FHC is fraction of ICE-COVERED area in this elevation class
    // Therefore, FHC=0 for the sea portion of the SeaLand Segment
    // NOT: fhc[sealand.base,_maskA] = 1.-fgice[_maskA]
    // We could do away with this EC altogether because it's not used.
    for (int iA_d=0; iA_d<dimA_global.dense_extent(); ++iA_d) {
        int iA_s = dimA_global.to_sparse(iA_d);

        fhcE2(sealand.base, iA_s) = 0.;
        undericeE2(sealand.base, iA_s) = 0;
        elevE2(sealand.base, iA_s) = 0.;
    };
// return;    // good
    // ihc=1: Ice portion of grid cell at mean for the ice portion
    // FHC is fraction of ICE-COVERED area in this elevation class
    // Therefore, FHC=1 for the land portion of the SeaLand Segment
    // NOT: fhc[sealand.base+1,_maskA] = fgice[_maskA]
    for (int iA_d=0; iA_d<dimA_global.dense_extent(); ++iA_d) {
        int iA_s = dimA_global.to_sparse(iA_d);

        fhcE2(sealand.base+1, iA_s) = sealand_mult * fgiceA(iA_s);
        undericeE2(sealand.base+1, iA_s) = UI_NOTHING;
        elevE2(sealand.base+1, iA_s) = elevA_global(iA_s);
    }

    // ---------- Segment 2: Full Elevation Classes
    for (auto ii=begin(AvE_global_e); ii != end(AvE_global_e); ++ii) {
        auto const iA = ii->index(0);
        auto const iE = ii->index(1);

        auto const iE_tuple(gcmA->indexingHC.index_to_tuple<int,2>(iE));
        auto const iA2 = iE_tuple[0];
        auto const ihc = iE_tuple[1];

        if (iA2 != iA) (*icebin_error)(-1,
            "Matrix is non-local: iA=%d, iE=%d, iA2=%d", (int)iA, (int)iE, (int)iA2);

        fhcE2(ec.base+ihc, iA) = ii->value() * ec_mult;
        undericeE2(ec.base+ihc, iA) = UI_ICEBIN;
    }

    for (int ihc=0; ihc<nhc_ice; ++ihc) {
    for (int iA=0; iA<nA; ++iA) {
        int ihcx = ec.base + ihc;

        if (ihcx < 0 || ihcx >= elevE2.extent(0)) (*icebin_error)(-1, "ihcx out of bounds: %d %d\n", ihcx, elevE2.extent(0));
        if (iA < 0 || iA >= elevE2.extent(1)) (*icebin_error)(-1, "ihcx out of bounds: %d %d\n", ihcx, elevE2.extent(1));

        elevE2(ihcx, iA) = gcmA->hcdefs[ihc];
    }}

    // ==================================================
    // 4) Fix bottom of atmosphere

    // Atmosphere sees mean of elevations over entire grid cell
    for (int iA_d=0; iA_d<dimA_global.dense_extent(); ++iA_d) {
        int iA_s = dimA_global.to_sparse(iA_d);

        zatmoA(iA_s) = elevE2(legacy.base, iA_s);
    }

}

/** This needs to be run at least once before matrices can be generated. */
void GCMCoupler_ModelE::update_topo(double time_s, bool initial_timestep)
{
    GCMRegridder_ModelE *gcmA = dynamic_cast<GCMRegridder_ModelE *>(&*gcm_regridder);

    // Transfer elemask and sigma parameters from each ice coupler
    std::vector<ElevMask<1>> elevmasks;
    std::vector<std::array<double,3>> sigmas;
    for (size_t sheet_index=0; sheet_index < gcmA->ice_regridders().index.size(); ++sheet_index) {
        IceCoupler *icec(&*ice_couplers[sheet_index]);
        elevmasks.push_back(ElevMask<1>(icec->elevmaskI, icec->maskI));
        sigmas.push_back(icec->sigma);
    }

    Topos &topoA(modele_inputs);
    icebin::modele::update_topo(
        gcmA, topoO_fname, elevmasks, sigmas,
        initial_timestep, gcm_params.hc_segments, gcm_params.primary_segment,
        topoA, foceanOm0);
}



// NetCDF Logging Stuff





}}
