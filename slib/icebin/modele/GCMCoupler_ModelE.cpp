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
#include <ibmisc/iostream.hpp>
#include <icebin/modele/GCMCoupler_ModelE.hpp>
#include <icebin/modele/grids.hpp>
#include <icebin/contracts/contracts.hpp>
#include <boost/filesystem.hpp>
#include <icebin/modele/GCMRegridder_ModelE.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/modele/topo.hpp>
#include <icebin/modele/merge_topo.hpp>
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
// Useful debugging printout functions; to be enabled if used.
#if 0
static void print_stuffE(GCMCoupler_ModelE const *self,
    VectorMultivec const &gcm_ivalsE_s,
    VarSet const &gcm_inputsE)
{
    printf("gcm_ivalsE_s.index.size() = %ld\n", gcm_ivalsE_s.index.size());
    printf("gcm_ivalsE_s.weights.size() = %ld\n", gcm_ivalsE_s.weights.size());
    int n=0;
    for (size_t i=0; i<gcm_ivalsE_s.size(); ++i) {
        long iE = gcm_ivalsE_s.index[i];
        auto ijk(self->gcm_regridder->indexing(GridAE::E).index_to_tuple<int,3>(iE));    // zero-based, alphabetical order
        if (ijk[0] == 0) {
            double const *vals = &gcm_ivalsE_s.vals[i*gcm_ivalsE_s.nvar];
            if (std::isnan(vals[1])) continue;
            printf("    %ld j=%d ihc=%d (w=%g):", iE, ijk[1], ijk[2], gcm_ivalsE_s.weights[i]);
            for (int j=0; j<gcm_ivalsE_s.nvar; ++j) printf(" %g", vals[j]);
            printf("\n");
            if (++n >= 10) return;
        }
    }
}

static void print_stuffE(GCMCoupler_ModelE const *self, GCMInput const &out)
{
    VectorMultivec const &gcm_ivalsE_s(out.gcm_ivalss_s[(int)IndexAE::ETOPO]);
    VarSet const &gcm_inputsE(self->gcm_inputs[(int)IndexAE::ETOPO]);
    print_stuffE(self, gcm_ivalsE_s, gcm_inputsE);
}



static void print_stuffA(GCMCoupler_ModelE const *self,
    VectorMultivec const &gcm_ivalsA_s,
    VarSet const &gcm_inputsA)
{
    printf("gcm_ivalsA_s.index.size() = %ld\n", gcm_ivalsA_s.index.size());
    printf("gcm_ivalsA_s.weights.size() = %ld\n", gcm_ivalsA_s.weights.size());
    int n=0;
    for (size_t i=0; i<gcm_ivalsA_s.size(); ++i) {
        long iA = gcm_ivalsA_s.index[i];
        auto ij(self->gcm_regridder->indexing(GridAE::A).index_to_tuple<int,2>(iA));    // zero-based, alphabetical order
        if (ij[0] == 0) {
            double const *vals = &gcm_ivalsA_s.vals[i*gcm_ivalsA_s.nvar];
            printf("    %ld j=%d (w=%g):", iA, ij[1], gcm_ivalsA_s.weights[i]);
            for (int j=0; j<gcm_ivalsA_s.nvar; ++j) printf(" %g", vals[j]);
            printf("\n");
            if (++n >= 10) return;
        }
    }
}

static void print_stuffA(GCMCoupler_ModelE const *self, GCMInput const &out)
{

    VectorMultivec const &gcm_ivalsA_s(out.gcm_ivalss_s[(int)IndexAE::ATOPO]);
    VarSet const &gcm_inputsA(self->gcm_inputs[(int)IndexAE::ATOPO]);
    print_stuffA(self, gcm_ivalsA_s, gcm_inputsA);
}
#endif
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

    scalars.add("by_dt", nan, "s-1", "", 0, "Inverse of coupling timestep");


    // Set up gcm_inputs array   (see IndexAE:  enum class { A, E, ATOPO, ETOPO, COUNT} IndexAE;)
    gcm_inputs_grid = std::vector<char>(indexae_grid);
    for (int i=0; i<(int)IndexAE::COUNT; ++i) {
        gcm_inputs.push_back(VarSet());
        gcm_ivalssA.push_back({});
        gcm_ivalssE.push_back({});
    }
}
// -----------------------------------------------------
static std::string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot); 
}

/** Returns the name of an ice sheet restart file, based on the (root
of) the GM restart file, i.e. sans .nc extension */
static std::string sheet_rsf(
    std::string const &gcm_rsf_root,
    std::string const &sheet_name)
{
    return strprintf("%s-%s.nc",
        gcm_rsf_root.c_str(), sheet_name.c_str());
}
// -----------------------------------------------------
IceCoupler::Params GCMCoupler_ModelE::make_ice_coupler_params(
    std::string const &sheet_name) const
{
    IceCoupler::Params params;

    // -----------------------
    // Determine name of the restart file used by ModelE
    // See MODELE.f for ISTART parameter meanings
    boost::filesystem::path gcm_rsf;
    switch (rdparams->istart) {
        case 2:    // observed start (cold start); no restart file
            params.rsf_fname = "";
            goto end_rsf;   // Done!
        case 8:
            gcm_rsf = "AIC";
            break;
        case 11:    // Restart from checkpoint file
            gcm_rsf = "fort.1.nc";
            break;
        case 12:    // Restart from checkpoint file
            gcm_rsf = "fort.2.nc";
            break;
        case 14:    // Restart from checkpoint file
            gcm_rsf = "fort.4.nc";
            break;
        default:
            (*icebin_error)(-1, "Unsupported ISTART value: %d", rdparams->istart);
    }

    // Follow the symlink
    // (expected for AIC and fort.4.nc)
    if (boost::filesystem::is_symlink(gcm_rsf))
        gcm_rsf = boost::filesystem::read_symlink(gcm_rsf);

    // Determine ice sheet restart filename
    params.rsf_fname = strprintf("%s-%s.nc",
        remove_extension(gcm_rsf.string()).c_str(),
        sheet_name.c_str());
end_rsf:
    // -----------------------

    return params;
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
    get_or_put_att(config_info, ncio_config.rw, "global_ec", global_ecO_fname);  // Must contain EvA at the very least

    /** EOpvAOp matrix for global (base) ice */
    ibmisc::ZArray<int,double,2> EOpvAOp_base;

    // Let's assume that the EOvAO matrix is included in topoO_fname

    // Replace the GCMRegridder with a wrapped version that understands
    // the ocean-vs-atmosphere grid complexity of ModelE
    gcm_regridder.reset(
        new modele::GCMRegridder_WrapE(    // allocates foceanOp and foceanOm
            std::unique_ptr<modele::GCMRegridder_ModelE>(
                new modele::GCMRegridder_ModelE(
                    global_ecO_fname,
                    std::dynamic_pointer_cast<GCMRegridder_Standard>(gcm_regridder)))));
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

    self->rdparams = &_rdparams;
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
printf("nhc_gcm2 %d\n", nhc_gcm);
    icebin_base_hc = 0;
    nhc_ice = self->gcm_regridder->nhc();
}

int GCMCoupler_ModelE::_read_nhc_gcm()
{
    // Get the name of the grid file
    std::string grid_fname;
    std::string global_ec_fname;
    {
        NcIO ncio_config(gcm_params.icebin_config_fname, NcFile::read);
        auto config_info(get_or_add_var(ncio_config, "m.info", "int", {}));
        get_or_put_att(config_info, ncio_config.rw, "grid", grid_fname);
        get_or_put_att(config_info, ncio_config.rw, "global_ec", global_ec_fname);
    }

    // Open the grid file to get the NHC

    // EC's for ice to be coupled to the ice sheet
    int nhc_grid;
    {NcIO ncio(grid_fname, 'r');
        nhc_grid = ncio.nc->getDim("m.nhc").getSize();
    }

    // EC's from global ice
    int nhc_global;
    {NcIO ncio(global_ec_fname, 'r');
        nhc_global = ncio.nc->getDim("nhc").getSize();
    }

    // Total ECs in the AvE matrix generated by merge_topo.cpp
    int const nhc_ice = nhc_grid + nhc_global;

    // Add another "EC" for the "land" segment in GCM runs.
    // This EC is NOT covered by any matrices.
    int const nhc_gcm = nhc_ice + 1;

    return nhc_gcm;

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
    std::string const name(name_f, name_len);
    self->gcm_constants.set(
        name, val,
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
char const *ncunits_f, int ncunits_len,
double mm, double bb,
char const *long_name_f, int long_name_len)
{
    std::string field_name(field_name_f, field_name_len);
    std::string units(units_f, units_len);
    std::string ncunits(ncunits_f, ncunits_len);
    std::string long_name(long_name_f, long_name_len);
    std::unique_ptr<blitz::Array<double,3>> var(
        new blitz::Array<double,3>(f_to_c(var_f.to_blitz())));

    unsigned int flags = 0;

    static double const xnan = std::numeric_limits<double>::quiet_NaN();
    self->gcm_outputsE.add(
        field_name, xnan, units, ncunits, mm, bb, flags, long_name);

    self->gcm_ovalsE.push_back(std::move(var));
}
// -----------------------------------------------------
/** @para var_nhc Number of elevation points for this variable.
 (equal to 1 for atmosphere variables, or nhc for elevation-grid variables)
@param return: Start of this variable in the gcm_inputs_local array (Fortran 1-based index) */
extern "C"
void gcmce_add_gcm_inputa(
GCMCoupler_ModelE *self,
int index_ae,
F90Array<double, 2> &var_f,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
char const *ncunits_f, int ncunits_len,
double mm, double bb,
bool initial,    // bool
char const *long_name_f, int long_name_len)
{
    if (self->gcm_inputs_grid[index_ae] != 'A') (*icebin_error)(-1,
        "gcmce_add_gcm_inputa() trying to add to VarSet on grid '%c'",
       self-> gcm_inputs_grid[index_ae]);

    std::string field_name(field_name_f, field_name_len);
    std::string units(units_f, units_len);
    std::string ncunits(ncunits_f, ncunits_len);
    std::string long_name(long_name_f, long_name_len);

    if (var_f.base == NULL) (*icebin_error)(-1,
        "gcmce_add_gcm_inputa() trying to add unallocated variable %s",
        field_name.c_str());
    std::unique_ptr<blitz::Array<double,2>> var(
        new blitz::Array<double,2>(f_to_c(var_f.to_blitz())));

    unsigned int flags = 0;
    if (initial) flags |= contracts::INITIAL;

    static double const xnan = std::numeric_limits<double>::quiet_NaN();
    self->gcm_inputs[index_ae].add(
        field_name, xnan, units, ncunits, mm, bb, flags, long_name);

    self->gcm_ivalssA[index_ae].push_back(std::move(var));
}
// -----------------------------------------------------
extern "C"
void gcmce_add_gcm_inpute(
GCMCoupler_ModelE *self,
int index_ae,
F90Array<double, 3> &var_f,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
char const *ncunits_f, int ncunits_len,
double mm, double bb,
int initial,    // bool
char const *long_name_f, int long_name_len)
{
    if (self->gcm_inputs_grid[index_ae] != 'E') (*icebin_error)(-1,
        "gcmce_add_gcm_inpute() trying to add to VarSet on grid '%c'",
        self->gcm_inputs_grid[index_ae]);

    std::string field_name(field_name_f, field_name_len);
    std::string units(units_f, units_len);
    std::string ncunits(ncunits_f, ncunits_len);
    std::string long_name(long_name_f, long_name_len);

    if (var_f.base == NULL) (*icebin_error)(-1,
        "gcmce_add_gcm_inpute() trying to add unallocated variable %s",
        field_name.c_str());
    std::unique_ptr<blitz::Array<double,3>> var(
        new blitz::Array<double,3>(f_to_c(var_f.to_blitz())));

    unsigned int flags = 0;
    if (initial) flags |= contracts::INITIAL;

    static double const xnan = std::numeric_limits<double>::quiet_NaN();
    self->gcm_inputs[index_ae].add(
        field_name, xnan, units, ncunits, mm, bb, flags, long_name);

    self->gcm_ivalssE[index_ae].push_back(std::move(var));
}
// -----------------------------------------------------
// ===========================================================
// Called from LISheetIceBin::io_rsf()   (warm starts)

extern "C"
void gcmce_write_rsf(GCMCoupler_ModelE *self,
    char *modele_fname_c, int modele_fname_n)
{
    // if (!self.am_i_root()) return;

    // modele_fname is filename of main ModelE restart file
    // We must modify it to produce name(s) of PISM restart file(s)
    // (one per ice sheet)
    std::string modele_fname(modele_fname_c, modele_fname_n);
    std::string modele_root(modele_fname);    // Already no .nc
    // std::string modele_root(remove_extension(modele_fname));

    // Save each PISM state...
    for (size_t sheetix=0; sheetix < self->ice_couplers.size(); ++sheetix) {
        auto &ice_coupler(self->ice_couplers[sheetix]);

        // Derive name for per-ice sheet PISM restart file
        ice_coupler->write_rsf(sheet_rsf(
            modele_root, ice_coupler->name()));
    }
}



// ===========================================================
// Called from LISheetIceBin::model_start()

extern "C"
void gcmce_model_start(GCMCoupler_ModelE *self, bool cold_start, int yeari, int itimei, double dtsrc)
{
    printf("BEGIN gcmce_model_start() cold_start=%d yeari=%d, itimei=%d, dtsrc=%g\n", cold_start, yeari, itimei, dtsrc);

    // This will cold-start initial conditions of the dynamic ice model
    self->dtsrc = dtsrc;

    // Call superclass model_start()
    double const time_s = itimei * dtsrc;
    self->model_start(cold_start,
        ibmisc::Datetime(yeari,1,1), time_s);

    // NOTE: Not needed because these things MUST be properly computed
    //       in the initial conditions TOPO file loaded by ModelE
    // b) Compute fhc, elevE
    // c) Compute ZATMO, FGICE, etc.
    //    self->update_topo(time_s);    // initial_timestep=true


    if (cold_start) {
        // d) Sync with dynamic ice model
        // This receives info back from ice model
        // (for warm start, the infor was already saved in a restart file)
        gcmce_couple_native(self, itimei, false,    // run_ice=false
            nullptr, nullptr, nullptr);    // !run_ice ==> no E1vE0c to return
    }

    printf("END gcmce_model_start()\n");
}


/** Helper function: splits a single GCMInput struct into per-domain GCMInput structs */
std::vector<GCMInput> split_by_domain(
    GCMInput const &out,
    DomainDecomposer_ModelE const &domainsA,
    DomainDecomposer_ModelE const &domainsE)
{
    using namespace spsparse;

    // Put domain decomposers in a nice array
    std::array<DomainDecomposer_ModelE const *, (int)IndexAE::COUNT> domainsAE
        {&domainsA, &domainsE, &domainsA, &domainsE};
    int ndomains = domainsA.size();    // Number of MPI domains

    // Construct the output
    std::vector<int> const &nvar(out.nvar());

    std::vector<GCMInput> outs;
    outs.reserve(ndomains);
    for (size_t i=0; i<ndomains; ++i)
        outs.push_back(GCMInput(nvar));

    size_t strideb = sizeof(GCMInput);

    // Split each element of parallel sparse vectors
    for (int iAE=0; iAE != (int)IndexAE::COUNT; ++iAE) {
        auto &gcm_ivalsX(out.gcm_ivalss_s[iAE]);
        DomainDecomposer_ModelE const &domainsX(*domainsAE[iAE]);

        for (size_t i=0; i<out.gcm_ivalss_s[iAE].index.size(); ++i) {
            // Convert from index to tuple
            long const iA(gcm_ivalsX.index[i]);

            // Figure out MPI domain of the resulting tuple
            int idomain = domainsX.get_domain(iA);

            // Add our values to the appropriate MPI domain
            outs[idomain].gcm_ivalss_s[iAE].add(
                iA, &gcm_ivalsX.vals[i*nvar[iAE]], gcm_ivalsX.weights[i]);
        }
    }

    // Split E1vE0
    // (E1vE0 is not set the first time around; in that case, shape = (-1,-1)
    if (out.E1vE0c.shape()[0] != -1) {
        auto shape(out.E1vE0c.shape());

        // Initialize output matrices
        for (size_t i=0; i<outs.size(); ++i) outs[i].E1vE0c.set_shape(shape);

        // Copy the main matrix
        // Works for matrix in A or E
        for (auto &tp : out.E1vE0c.tuples) {
            int const domain = domainsE.get_domain(tp.index(0));
            outs[domain].E1vE0c.add(tp.index(), tp.value());
        }
    }
    return outs;
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
bool run_ice,    // if false, only initialize
// https://stackoverflow.com/questions/30152073/how-to-pass-c-pointer-to-fortran
// Return the E1vE0 matrix here
int **E1vE0c_indices_p,
double **E1vE0c_values_p,
int *E1vE0c_nele)
{
    double time_s = itime * self->dtsrc;

    // Fill it in...
    VectorMultivec gcm_ovalsE_s(self->gcm_outputsE.size());
    std::vector<double> val(self->gcm_outputsE.size());    // Temporary

    auto const &indexingA(self->gcm_regridder->agridA->indexing);
    auto const &indexingE(self->gcm_regridder->indexingE);

    // domain uses alphabetical order, 0-based indexing...
    const auto base_hc(self->gcm_params.icebin_base_hc);
    const auto nhc_ice(self->gcm_regridder->nhc());

    // Lookup underice (C-style order and indexing base)
    int const i_underice = self->gcm_inputs[(int)IndexAE::ETOPO].index.at("underice_d");
    blitz::Array<double,3> &underice(*self->gcm_ivalssE[(int)IndexAE::ETOPO][i_underice]);
    

    // Get values to send to IceBin
    auto &domainA(self->domainA);
printf("domainA size=%ld base_hc=%d  nhc_ice=%d\n", domainA.data.size(), base_hc, nhc_ice);

    for (int ihc=base_hc; ihc < base_hc + nhc_ice; ++ihc) {
        const int ihc_ice = ihc - base_hc;   // Use only IceBin HC's, zero-based indexing
        if (ihc_ice < 0) (*icebin_error)(-1,
            "ihc_ice cannot be <0: %d = %d - %d - 1\n", ihc_ice, ihc, base_hc);

        // Iterate over just this MPI rank's chunk
        for (int j=domainA[1].begin; j < domainA[1].end; ++j) {
        for (int i=domainA[0].begin; i < domainA[0].end; ++i) {
            // i,j are 0-based indexes.
            if (underice(ihc,j,i) == UI_LOCALICE || underice(ihc,j,i) == UI_GLOBALICE) {

                long iE_s = indexingE.tuple_to_index(
                    make_array(i, j, ihc_ice));

                if (iE_s < 0) (*ibmisc_error)(-1,
                    "iE_s=%ld (from %d %d %d), it should not be negative\n", iE_s, i, j, ihc_ice);

                for (unsigned int ivar=0; ivar<self->gcm_outputsE.size(); ++ivar) {
                    val[ivar] = (*self->gcm_ovalsE[ivar])(ihc,j,i);
                }

                gcm_ovalsE_s.add({iE_s}, &val[0], 1.0);
            }    // if UI_LOCALICE or UI_GLOBALICE
        }}
    }


    // Gather it to root
    // boost::mpi::communicator &gcm_world(world);
    // Init our output struct based on number of A and E variables.
    std::vector<int> sizes;
    for (VarSet const &vs : self->gcm_inputs) sizes.push_back(vs.size());
    GCMInput out(sizes);

    if (self->am_i_root()) {
        // =================== MPI ROOT =============================
        std::vector<VectorMultivec> every_gcm_ovalsE_s;
        boost::mpi::gather(self->gcm_params.world, gcm_ovalsE_s,
            every_gcm_ovalsE_s, self->gcm_params.gcm_root);

        // Concatenate coupler inputs
        VectorMultivec gcm_ovalsE_s(concatenate(every_gcm_ovalsE_s));

        // Couple on root!
        // out contains GLOBAL output for all MPI ranks
        out = self->couple(time_s, gcm_ovalsE_s, run_ice);  // move semantics

        // Split up the output (and 
        std::vector<GCMInput> every_outs(
            split_by_domain(out, *self->domains, *self->domains));

        // Scatter!
        boost::mpi::scatter(self->gcm_params.world, every_outs, out, self->gcm_params.gcm_root);


    } else {
        // =================== NOT MPI ROOT =============================
        // Send our input to root
        boost::mpi::gather(self->gcm_params.world, gcm_ovalsE_s, self->gcm_params.gcm_root);

        // Let root do the work...
        // update_topo() is built into this
        self->couple(time_s, gcm_ovalsE_s, run_ice);

        // Receive our output back from root
        boost::mpi::scatter(self->gcm_params.world, out, self->gcm_params.gcm_root);
    }

    // 1. Copies values back into modele.gcm_ivals from scatterd MPI stuff
    self->apply_gcm_ivals(out);

    auto const &indexingHC(self->gcm_regridder->indexingHC);
    // Copy E1vE0 matrix back to Fortran
    self->E1vE0c.clear();
    if (run_ice) {
        for (auto ii(out.E1vE0c.begin()); ii != out.E1vE0c.end(); ++ii) {
            long const iE1(ii->index(0));
            long const iE0(ii->index(1));
            
            auto ijk0(indexingHC.index_to_tuple<int,2>(iE0));
            auto ijk1(indexingHC.index_to_tuple<int,2>(iE1));
            if (ijk0[0] != ijk1[0]) (*icebin_error)(-1,
                "ijk0 (%d) and ijk1 (%d) must match (%g)",
                ijk0[0], ijk1[0], ii->value());
            auto ij(indexingA.index_to_tuple<int,2>(ijk0[0]));

            int const ihc0 = ijk0[1];
            int const ihc1 = ijk1[1];

            // +1: Convert to Fortran 1-based Indexing
            self->E1vE0c.indices.push_back(ij[0]+1);
            self->E1vE0c.indices.push_back(ij[1]+1);
            self->E1vE0c.indices.push_back(ihc1+1);
            self->E1vE0c.indices.push_back(ihc0+1);
            self->E1vE0c.values.push_back(ii->value());
        }

        // Send E1vE0 back to ModelE in Fortran
        // See https://stackoverflow.com/questions/30152073/how-to-pass-c-pointer-to-fortran
        // self->E1vE0c.export_matrix(*E1vE0c_indices_p, *E1vE0c_values_p, *E1vE0c_nele);
        *E1vE0c_nele = self->E1vE0c.values.size();
        *E1vE0c_indices_p = &self->E1vE0c.indices[0];
        *E1vE0c_values_p = &self->E1vE0c.values[0];


    }
}
// =======================================================

/** Called from MPI rank.  Copies output of coupling back into
appropriate dense-indexing ModelE variables. */
void GCMCoupler_ModelE::apply_gcm_ivals(GCMInput const &out)
{
    printf("BEGIN GCMCoupler_ModelE::apply_gcm_ivals\n");
    auto nvar(out.nvar());    // A and E
    const auto nhc_ice(this->gcm_regridder->nhc());

    // Write to here...
    for (int iAE=0; iAE<GridAE::count; ++iAE) {
        if (nvar[iAE] < 0) (*icebin_error)(-1,
            "nvar[%d]=%d < 0, it should not be\n", iAE, nvar[iAE]);
    }

    for (size_t index_ae=0; index_ae < gcm_inputs.size(); ++index_ae) {
    switch(gcm_inputs_grid[index_ae]) {
        case 'A' : {
            VectorMultivec const &gcm_ivalsA_s(out.gcm_ivalss_s[index_ae]);    // src
            std::vector<std::unique_ptr<blitz::Array<double,2>>> &gcm_ivalsA(
                gcm_ivalssA[index_ae]); // dest

            // ================================= A

            // Update gcm_ivalA variables...
            // Read from here...
            if (gcm_ivalsA.size() != gcm_inputs[index_ae].size()) (*icebin_error)(-1,
                "gcm_ivalsA is wrong size: %ld vs. %ld", gcm_ivalsA.size(), gcm_inputs[index_ae].size());

            // Clear output; but only for gridcells we want to touch
            int const nvar = gcm_inputs[index_ae].size();
            for (size_t ix=0; ix<gcm_ivalsA_s.size(); ++ix) {    // Iterate through elements of parallel arrays
                long const iA = gcm_ivalsA_s.index[ix];
                auto ij(gcm_regridder->indexing(GridAE::A).index_to_tuple<int,2>(iA));    // zero-based, alphabetical order
                int const i = ij[0];
                int const j = ij[1];
                for (int ivar=0; ivar<nvar; ++ivar) {
                    (*gcm_ivalsA[ivar])(j,i) = 0;
                }
            }

            // Create (inverse of) summed weights
            blitz::Array<double,1> sA_s(gcm_regridder->nA());
            gcm_ivalsA_s.to_dense_scale(sA_s);

            // Copy sparse arrays to output
            for (size_t ix=0; ix<gcm_ivalsA_s.size(); ++ix) {    // Iterate through elements of parallel arrays
                long const iA = gcm_ivalsA_s.index[ix];
                auto ij(gcm_regridder->indexing(GridAE::A).index_to_tuple<int,2>(iA));    // zero-based, alphabetical order
                int const i = ij[0];
                int const j = ij[1];
                for (int ivar=0; ivar<nvar; ++ivar) {
                    // Original Fortran (partial) arrays have been converted to
                    // C++ order and 0-based indexing; see gcmce_add_gcm_inputa()
                    (*gcm_ivalsA[ivar])(j,i) +=
                        gcm_ivalsA_s.vals[ix*nvar + ivar] * sA_s(iA);
                }
            }
        } break;
        case 'E' : {
            // ================================= E
            VectorMultivec const &gcm_ivalsE_s(out.gcm_ivalss_s[index_ae]);    // src
            std::vector<std::unique_ptr<blitz::Array<double,3>>>  &gcm_ivalsE(
                gcm_ivalssE[index_ae]); // dest

            // Update gcm_ivalE variables...
            // Read from here...
            if (gcm_ivalsE.size() != gcm_inputs[index_ae].size()) (*icebin_error)(-1,
                "gcm_ivalsE is wrong size: %ld vs. %ld (index_ae = %d)", gcm_ivalsE.size(), gcm_inputs[index_ae].size(), index_ae);

            // Clear output: because non-present elements in sparse gcm_ivalsA are 0
            int const nvar = gcm_inputs[index_ae].size();
            for (size_t ix=0; ix<gcm_ivalsE_s.size(); ++ix) {
                long const iE = gcm_ivalsE_s.index[ix];
                auto ijk(gcm_regridder->indexing(GridAE::E).index_to_tuple<int,3>(iE));
                int const i = ijk[0];
                int const j = ijk[1];
                int const ihc_ice = ijk[2];    // zero-based, just EC's known by ice model
                int const ihc_gcm = ihc_ice;

                for (int ivar=0; ivar<nvar; ++ivar) {
                    // Original Fortran (partial) arrays have been converted to
                    // C++ order and 0-based indexing; see gcmce_add_gcm_inputa()
                    (*gcm_ivalsE[ivar])(blitz::Range(0,nhc_ice-1),j,i) = 0;
                }
            }

            // Create (inverse of) summed weights
            blitz::Array<double,1> sE_s(gcm_regridder->nE());
            gcm_ivalsE_s.to_dense_scale(sE_s);

            // Copy sparse arrays to output
            for (size_t ix=0; ix<gcm_ivalsE_s.size(); ++ix) {
                long const iE = gcm_ivalsE_s.index[ix];
                auto ijk(gcm_regridder->indexing(GridAE::E).index_to_tuple<int,3>(iE));
                int const i = ijk[0];
                int const j = ijk[1];
                int const ihc_ice = ijk[2];    // zero-based, just EC's known by ice model
                int const ihc_gcm = ihc_ice;

                for (int ivar=0; ivar<nvar; ++ivar) {
                    // Original Fortran (partial) arrays have been converted to
                    // C++ order and 0-based indexing; see gcmce_add_gcm_inputa()
                    (*gcm_ivalsE[ivar])(ihc_gcm,j,i) +=
                        gcm_ivalsE_s.vals[ix*gcm_inputs[index_ae].size() + ivar] * sE_s(iE);
                }
            }
        } break;
        default :
            (*icebin_error)(-1, "Illegal gcm_inputs_grid of '%c'", gcm_inputs_grid[index_ae]);
    }}

    printf("END GCMCoupler_ModelE::apply_gcm_ivals\n");
}
// ============================================================================
// Update TOPO file during a coupled run

#if 0
/** Produce regridding matrices for this setup.
(To be run on root MPI node)
Needs to include _foceanAOp and _foceanAOm */
std::unique_ptr<RegridMatrices_Dynamic> GCMCoupler_ModelE::regrid_matrices(    // virtual
    int sheet_index,
    blitz::Array<double,1> const &elevmaskI,
    RegridParams const &params) const
{
    GCMRegridder_ModelE const *gcmA(
        dynamic_cast<GCMRegridder_WrapE *>(&*gcm_regridder)->gcmA.get());

    return gcmA->regrid_matrices(
        sheet_index,
        _foceanAOp, _foceanAOm,
        elevmaskI, params);
}
#endif


void GCMCoupler_ModelE::update_topo(
double time_s,    // Simulation time
bool run_ice,     // false for initialization
std::vector<blitz::Array<double,1>> const &emI_lands,
std::vector<blitz::Array<double,1>> const &emI_ices,
// ---------- Input & Output
// Write: gcm_ivalss_s[IndexAE::ATOPO], gcm_ivalss_s[IndexAE::ETOPO]
GCMInput &out,
TupleListLT<1> &wEAm_base)   // Clear; then store wEAm in here
{
    auto const &indexingA(gcm_regridder->agridA->indexing);
    auto const &indexingE(gcm_regridder->indexingE);

    GCMRegridder_WrapE *gcmW(
        dynamic_cast<GCMRegridder_WrapE *>(&*gcm_regridder));
    GCMRegridder_ModelE const *gcmA(gcmW->gcmA.get());

    // Read from TOPOO file
    ibmisc::ArrayBundle<double,2> topoo(topoo_bundle(BundleOType::MERGEO, topoO_fname));
        auto &foceanOp(topoo.array("FOCEANF"));
        auto &fgiceOp(topoo.array("FGICEF"));
        auto &zatmoOp(topoo.array("ZATMOF"));
        auto &foceanOm(topoo.array("FOCEAN"));
        auto &flakeOm(topoo.array("FLAKE"));
        auto &fgrndOm(topoo.array("FGRND"));
        auto &fgiceOm(topoo.array("FGICE"));
        auto &zatmoOm(topoo.array("ZATMO"));
        auto &zlakeOm(topoo.array("ZLAKE"));
        auto &zicetopO(topoo.array("ZICETOP"));
        auto &zland_minO(topoo.array("ZLAND_MIN"));
        auto &zland_maxO(topoo.array("ZLAND_MAX"));
        blitz::Array<int16_t,2> mergemaskO(zicetopO.extent());

    // Allocate space for TOPOA output variables (variables in TOPOA file)
    TopoABundles topoa(topoo, gcmA->hspecA(), nhc_gcm());
        auto &foceanA(topoa.a.array("focean"));
        auto &flakeA(topoa.a.array("flake"));
        auto &fgrndA(topoa.a.array("fgrnd"));
        auto &fgiceA(topoa.a.array("fgice"));
        auto &zatmoA(topoa.a.array("zatmo"));
        auto &hlakeA(topoa.a.array("hlake"));
        auto &zicetopA(topoa.a.array("zicetop"));
        auto &zland_minA(topoa.a.array("zland_min"));
        auto &zland_maxA(topoa.a.array("zland_max"));
        auto &mergemaskA(topoa.a_i.array("mergemask"));
        if (!run_ice) mergemaskA0.reference(mergemaskA);  // no mergemask on initialization

        auto &fhc(topoa.a3.array("fhc"));
        auto &elevE(topoa.a3.array("elevE"));

        auto &underice_i(topoa.a3_i.array("underice"));

    std::vector<std::string> errors;

    // ------------------ Create merged TOPOO file (in RAM)
    // We need correctA=true here to get FOCEANF, etc.
    // merge_topoO() merges PISM ice sheets (emI_ices, etc) into foceanOp / foceanOm
    merge_topoO(
        foceanOp, fgiceOp, zatmoOp,  // Our own bundle
        foceanOm, flakeOm, fgrndOm, fgiceOm, zatmoOm, zicetopO,
        zland_minO, zland_maxO, mergemaskO, &*gcmA->gcmO,
        RegridParams(false, true, {0.,0.,0.}),  // (scale, correctA, sigma)
        emI_lands, emI_ices, gcmA->specO().eq_rad, errors);

    // Copy FOCEAN to internal GCMRegridder_WrapE state
    gcmW->foceanOp = reshape1(foceanOp);
    // Only copy foceanAOm the first time.  It NEVER changes after that.
    if (!run_ice) {
        gcmW->foceanOm = reshape1(foceanOm);
    }

    // Print sanity check errors to STDERR
    for (std::string const &err : errors) fprintf(stderr, "ERROR: %s\n", err.c_str());

    if (errors.size() > 0) (*icebin_error)(-1,
        "Errors in TOPO merging or regridding; halting!");

    // ---------------- Compute AAmvEAm (requires merged foceanOp, foceanOm)
    long offsetE;  // Offset (in sparse E space) added to base EC indices
    linear::Weighted_Tuple AAmvEAm(gcmA->global_AvE(
        emI_lands, emI_ices, reshape1(foceanOp), reshape1(foceanOm),
        true, offsetE));    // scale=true

    // ---------------- Compute wAEm_base (weight of JUST base-ice ECs)
    // Base ice ECs are distinguished because they've been offsetted ("stacked")
    // in compute_EOpvAOp_merged() (merge_topo.cpp; squash_ecs=false)
    wEAm_base.clear();
    for (auto ii=AAmvEAm.Mw.begin(); ii != AAmvEAm.Mw.end(); ++ii) {
        if (ii->index(0) >= offsetE) {  // Only copy weight of base ECs
            wEAm_base.add(ii->index(), ii->value());
        }
    }

    // ---------------- Create TOPOA file (in RAM)
    std::vector<int16_t> underice_hc;
    for (int ihc=0; ihc<gcmA->nhc(); ++ihc) underice_hc.push_back(gcmA->underice(ihc));

    std::vector<std::string> errors2(make_topoA(
        foceanOm, flakeOm, fgrndOm, fgiceOm, zatmoOm, zlakeOm, zicetopO,
        zland_minO, zland_maxO, mergemaskO,
        gcmA->hspecO(), gcmA->hspecA(), gcmA->indexingHC, gcmA->hcdefs(), underice_hc,
        AAmvEAm,
        foceanA, flakeA, fgrndA, fgiceA, zatmoA, hlakeA, zicetopA,
        zland_minA, zland_maxA, mergemaskA,
        fhc, elevE, underice_i));

    // Print sanity check errors to STDERR
    for (std::string const &err : errors2) fprintf(stderr, "ERROR: %s\n", err.c_str());

    if (errors2.size() > 0) (*icebin_error)(-1,
        "Errors in TOPO merging or regridding; halting!");

    // ------------------ Deallocate stuff we no longer need

    // Compute wEAm_base
    AAmvEAm.clear();
    topoo.free();

    // ------------------- Convert underice to double and store in topoa
    auto &hspecA(gcmA->hspecA());
    auto const nhc(gcmA->nhc());
    std::array<int,3> shape3 {nhc, hspecA.jm, hspecA.im};
    topoa.a3.add("underice_d", {});
    topoa.a3.at("underice_d").allocate(true, shape3);
    auto &underice(topoa.a3.array("underice_d"));
    underice = underice_i;

    // ------------------ Pack TOPOA stuff into output variables
    {
        // See LISheetIceBin.F90 (modelE); same as variables found in classic TOPO file
        // These variables should be: "focean", "flake", "fgrnd", "fgice", "zatmo", "hlake"
        VarSet &gcm_inputsA(this->gcm_inputs[(int)IndexAE::ATOPO]);
        VectorMultivec &gcm_ivalsA_s(out.gcm_ivalss_s[(int)IndexAE::ATOPO]);

        // Map between contract and TOPOA variables
        std::vector<blitz::Array<double,2> *> iarrays;
        for (size_t k=0; k<gcm_inputsA.index.size(); ++k) {
            // k = index in contract
            // name = name of the contract variable
            auto const &name(gcm_inputsA[k].name);

            // iarrays_ix = Index in topoa.a
            size_t const iarrays_ix = topoa.a.index.at(name);
            iarrays.push_back(&topoa.a.array(iarrays_ix));
        }

        // Copy every element
        std::vector<double> val(gcm_ivalsA_s.nvar);
        for (int j=0; j<hspecA.jm; ++j) {
        for (int i=0; i<hspecA.im; ++i) {
            if (mergemaskA(j,i) == 0 && mergemaskA0(j,i) == 0) continue;

            for (int k=0; k<gcm_ivalsA_s.nvar; ++k) {
                val[k] = (*iarrays[k])(j,i);
            }
            auto ij = indexingA.tuple_to_index(std::array<int,2>{i,j});
            gcm_ivalsA_s.add(ij, val, 1.0);
        }}
    }

    // ------------------ Pack TOPOE stuff into output variables
    {
        // See LISheetIceBin.F90 (modelE); same as variables found in classic TOPO file
        // These variables should be: fhc, elevE, underice
        VarSet &gcm_inputsE(this->gcm_inputs[(int)IndexAE::ETOPO]);
        VectorMultivec &gcm_ivalsE_s(out.gcm_ivalss_s[(int)IndexAE::ETOPO]);

        // Map between contract and TOPOE variables
        std::vector<blitz::Array<double,3> *> iarrays;
        for (size_t k=0; k<gcm_inputsE.index.size(); ++k) {
            // k = index in contract
            // name = name of the contract variable
            auto const &name(gcm_inputsE[k].name);

            // iarrays_ix = Index in topoa.a
            int const iarrays_ix = topoa.a3.index.at(name);
            iarrays.push_back(&topoa.a3.array(iarrays_ix));
        }

        // Copy every element
        std::vector<double> val(gcm_ivalsE_s.nvar);
        for (int ihc=0; ihc<nhc; ++ihc) {
        for (int j=0; j<hspecA.jm; ++j) {
        for (int i=0; i<hspecA.im; ++i) {
            if (mergemaskA(j,i) == 0 && mergemaskA0(j,i) == 0) continue;

            for (int k=0; k<gcm_ivalsE_s.nvar; ++k) {
                val[k] = (*iarrays[k])(ihc,j,i);
            }
            auto ij = indexingE.tuple_to_index(std::array<int,3>{i,j,ihc});
            gcm_ivalsE_s.add(ij, val, 1.0);
        }}}
    }

    // Store this timestep's mergemask for next coupling time
    mergemaskA0.reference(mergemaskA);
}

// ----------------------------------------------------------------------
GCMInput GCMCoupler_ModelE::couple(
double time_s,        // Simulation time [s]
VectorMultivec const &gcm_ovalsE,
bool run_ice)    // if false, only initialize
{
    // Call superclass coupling for starters
    GCMInput out(this->GCMCoupler::couple(time_s, gcm_ovalsE, run_ice));

    // Nothing more to do unless we're root
    if (!gcm_params.am_i_root()) return out;

    // Run update_topo()
    TupleListLT<1> wEAm_base;  // set by update_topo()
    {
        std::vector<blitz::Array<double,1>> emI_lands, emI_ices;
        emI_ices.reserve(ice_couplers.size());
        emI_lands.reserve(ice_couplers.size());
        for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
            auto &ice_coupler(ice_couplers[sheetix]);

            emI_ices.push_back(ice_coupler->emI_ice);
            emI_lands.push_back(ice_coupler->emI_land);
        }

        update_topo(time_s, run_ice, emI_lands, emI_ices, out, wEAm_base);
    }


    // Log the results
    if (gcm_params.icebin_logging) {
        std::string fname = "gcm-in-" + sdate(time_s) + ".nc";
        NcIO ncio(fname, 'w');
        auto one_dims(get_or_add_dims(ncio, {"one"}, {1}));
        NcVar info_var = get_or_add_var(ncio, "info", ibmisc::get_nc_type<double>(), one_dims);
        info_var.putAtt("notes", "Elevation classes (HC) are just those known to IceBin.  No legacy or sea-land elevation classes included.");
        ncio_gcm_input(ncio, out, timespan, time_unit, "");
        ncio();
    }

    // ---------- Apply scaling to gcm_ivalsA_s, originally set in gcmce_add_xxx()

    // This converts (for example) ZATMO [m] (as needs to go into TOPO
    //     file) to ZATMO [m^2 s-2] (as ModelE wants to see internally)
    // It is more trouble-free to do this here, rather than in
    //     IceCoupler::couple() and update_topo()
    for (size_t index_ae=0; index_ae < gcm_inputs.size(); ++index_ae) {
        VarSet const &gcm_inputsA(gcm_inputs[index_ae]);
        int const nvar = gcm_inputsA.size();
        VectorMultivec &gcm_ivalsA_s(out.gcm_ivalss_s[index_ae]);

        for (size_t ix=0; ix<gcm_ivalsA_s.size(); ++ix) {    // Iterate through elements of parallel arrays
            for (int ivar=0; ivar<nvar; ++ivar) {
                double &val(gcm_ivalsA_s.vals[ix*nvar + ivar]);
                VarMeta const &gcm_input(gcm_inputsA.data[ivar]);
                val = val * gcm_input.mm + gcm_input.bb;
            }
        }
    }



printf("END GCMCoupler::couple()\n");
    return out;
}
// ------------------------------------------------------------


//TODO: Take care of FLAND, which was not computed by make_topoa


}}
