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

#pragma once

#include <boost/mpi.hpp>
#include <ibmisc/f90blitz.hpp>
#include <icebin/GCMCoupler.hpp>
#include <icebin/modele/GCMRegridder_ModelE.hpp>

namespace icebin {
namespace modele {

#if 0
// Don't need this for now...
BOOST_ENUM_VALUES( ModelE_CouplingType, int,
    /** GCM reports top T boundary condition to ice sheet.  This is
    always available. */
    (DIRICHLET_BC) (0)

    /** GCM reports energy fluxes at top of ice sheet.  This is only
    available on some ice models. */
    (NEUMANN_BC) (1)
);

#endif


// ---------------------------------------------
// Parameters read from the ModelE rundeck
// These have a peer in api_f.f90
static int const MAX_CHAR_LEN = 128;    // From Dictionary_mod.F90
struct ModelEParams
{
//    char icebin_segments[MAX_CHAR_LEN];
//    char ice_coupler_type[MAX_CHAR_LEN];    // DISMAL,PISM
//    double dtsrc;
    int dummy;    // Avoid zero-size struct
};
// ---------------------------------------------


class DomainDecomposer_ModelE {
    ibmisc::Domain domainA_global;
    size_t ndomain;
    blitz::Array<int,1> rank_of_j;    // indexing base=1
public:

    DomainDecomposer_ModelE(std::vector<int> const &endj, ibmisc::Domain const &_domainA_global);

    /** Number of domains */
    size_t size() const { return ndomain; }

    /** Returns the MPI rank of grid cell.  Works if ix is iA (atmosphere grid) or iE (elevation grid) */
    int get_domain(long ix) const {    // zero-based
        auto im_world(domainA_global[0].end);
        auto jm_world(domainA_global[1].end);

        int j = (ix / im_world) % jm_world;    // +0 for 0-based indexing
        return rank_of_j(j);
    }
};

#if 0
struct GCMInput_ModelE : public GCMInput
{

    enum class A {FOCEAN, FLAKE, FGRND, FGICE, ZATMO, ZLAKE, ZICETOP, NUM};
    enum class B {FHC, ELEV, UNDERICE, NUM};

    // TOPO file stuff
    //    * Computed by GCMCoupler_ModelE::couple()
    //    * SCALED: ready-to-use values
    std::array<VectorMultivec, GridAE::count> topoAE_s;
    std::array<std::vector<double>, GridAE::count> topoAE_weight_s;

    GCMInput_ModelE(std::array<int, GridAE::count> const &nvar) :
        GCMInput(nvar),   // From the contract
        topoAE_s({
            VectorMultivec(A::NUM),   // Add weight as last column
            VectorMultivec(E::NUM),
        })
    {}
};
#endif


class GCMCoupler_ModelE : public GCMCoupler
{
public:
    double dtsrc;
    ModelEParams rdparams;    // Params straight from the rundeck (came during init)

    /** On root: separate global stuff back into individual domains.
    Works for A and E grids. */
    std::unique_ptr<DomainDecomposer_ModelE> domains;

    // ================== ModelE Outputs
    // gcm_ovalsE[ovar](i, j, ihc)    Fortran-order 1-based indexing
    std::vector<std::unique_ptr<blitz::Array<double,3>>> gcm_ovalsE;

    // ================== ModelE Inputs
    // Variables borrowed from ModelE, used to return data to it.
    // All these variables are Fortran-order, 1-based indexing

    // References back to original ModelE-supplied MPI arrays
    // gcm_ivalsA[set_index][ivar](i, j, ihc)    C++-order 0 based ordering
    // Indexed the same as gcm_inputs.
    std::vector<std::vector<std::unique_ptr<blitz::Array<double,2>>>> gcm_ivalssA;
    std::vector<std::vector<std::unique_ptr<blitz::Array<double,3>>>> gcm_ivalssE;

    // Low and high indices for this MPI rank.
    // Indices are in Fortran order (im, jm) with zero-based indexing
    ibmisc::Domain domainA;
    // Low and high indices for global domain (Fortran order, 0-based)
    ibmisc::Domain domainA_global;

    /** Name of the Ocean-level TOPO file (output of modified Gary's
    program, sans ice sheets) */
    std::string topoO_fname;

    /** Name of file on ocean grid containing the EvA matrix for global (non-IceBin) ice. */
    std::string global_ecO_fname;

    /** ModelE ocean cover, on the Ocean grid, as seen by the ice
    model (sparse indexing).  Ocean grid cells can contain fractional
    ocean cover.  foceanAOp can change over the course of a ModelE
    run, as the ice model's ice extent changes.
    @see ModelE FOCEAN or FOCEN */
    blitz::Array<double,1> _foceanAOp;

    /** ModelE ocean cover, on the Ocean grid, as seen by ModelE
    (sparse indexing).  Cells are either all ocean (==1.0) or all
    continent (==0.0).  foceanAOm does NOT change over the course of a
    ModelE run, because the ModelE ocean is not able to change shape
    mid-run. */
    blitz::Array<double,1> _foceanAOm;


    /** Merge mask from last timestep */
    blitz::Array<int16_t,2> mergemaskA0;

public:
    virtual ~GCMCoupler_ModelE() {}

    // Called from LISnow::allocate()
    GCMCoupler_ModelE(GCMParams &&_params);

    /** @param am_i_root
        Call with true if calling from MPI root; false otherwise.
        The core coupling/regridding computation only runs on root.
        But other MPI ranks need to go along for the ride, assuming that
        the ice model uses MPI. */
    GCMInput couple(
        double time_s,        // Simulation time [s]
        VectorMultivec const &gcm_ovalsE,
        bool run_ice);    // if false, only initialize

    void _ncread(    // virtual
        ibmisc::NcIO &ncio_config,
        std::string const &vname);        // comes from this->gcm_params

    std::string locate_input_file(   // virtual
        std::string const &sheet_name,        // eg: greenland
        std::string const &file_name);        // eg: pism_Greenland_5km_v1.1.nc

    /**
    @param out Add TOPO stuff here (see IndexAE::TOPOA, IndexAE::TOPOE) */
    void update_topo(
    double time_s,    // Simulation time
    bool run_ice,     // false for initialization
    std::vector<blitz::Array<double,1>> const &emI_lands,
    std::vector<blitz::Array<double,1>> const &emI_ices,
    // ---------- Input & Output
    // Write: gcm_ivalss_s[IndexAE::ATOPO], gcm_ivalss_s[IndexAE::ETOPO]
    GCMInput &out,
    TupleListLT<1> &wEAm_base);   // Clear; then store wEAm in here

    int _read_nhc_gcm();

    /** Copies GCM inputs back to original GCM-supplied sparse input arrays */
    void apply_gcm_ivals(GCMInput const &out);

    // The gcmce_xxx() functions do not need to be declared here
    // because everything in this class is public.


//    // 1. Copies values back into modele_inputs.gcm_ivals
//    void update_gcm_ivals(GCMInput const &out);

};    // class GCMCouler_ModelE

// ===============================================================
// The "gcmce_*" interface used by Fortran ModelE

extern "C"
GCMCoupler_ModelE *gcmce_new(
    ModelEParams const &_rdparams,

    // Info about the global grid
    int im, int jm,

    // Info about the local grid (1-based indexing, closed ranges)
    int i0_f, int i1_f, int j0_f, int j1_f,

    // MPI Stuff
    MPI_Fint comm_f, int root);

/* Tells ModelE how many elevation classes it needs **/
extern "C"
void gcmce_hc_params(GCMCoupler_ModelE *self, int &nhc_gcm, int &icebin_base_hc, int &nhc_ice);

/** Set a single constant value in Icebin.  This is a callback, to be called
from ModelE's (Fortran code) constant_set::set_all_constants() */
/** Called from within LISheetIceBin::allocate() */
extern "C"
void gcmce_set_constant(
    GCMCoupler_ModelE *self,
    char const *name_f, int name_len,
    double val,
    char const *units_f, int units_len,
    char const *description_f, int description_len);

extern "C"
void gcmce_add_gcm_outpute(
GCMCoupler_ModelE *self,
ibmisc::F90Array<double, 3> &var_f,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
char const *long_name_f, int long_name_len);

/** @para var_nhc Number of elevation points for this variable.
 (equal to 1 for atmosphere variables, or nhc for elevation-grid variables)
@param return: Start of this variable in the gcm_inputs_local array (Fortran 1-based index) */
extern "C"
void gcmce_add_gcm_inputa(
GCMCoupler_ModelE *self,
int index_ae,
ibmisc::F90Array<double, 2> &var_f,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
bool initial,    // bool
char const *long_name_f, int long_name_len);

extern "C"
void gcmce_add_gcm_inpute(
GCMCoupler_ModelE *self,
int index_ae,
ibmisc::F90Array<double, 3> &var_f,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
int initial,    // bool
char const *long_name_f, int long_name_len);

extern "C"
void gcmce_reference_globals(
    GCMCoupler_ModelE *self,
    ibmisc::F90Array<double, 3> fhc,
    ibmisc::F90Array<int, 3> underice,
    ibmisc::F90Array<double, 3> elevE,
    ibmisc::F90Array<double, 2> focean,
    ibmisc::F90Array<double, 2> flake,
    ibmisc::F90Array<double, 2> fgrnd,
    ibmisc::F90Array<double, 2> fgice,
    ibmisc::F90Array<double, 2> zatmo);

extern "C"
void gcmce_io_rsf(GCMCoupler_ModelE *self,
    char *fname_c, int fname_n);

extern "C"
void gcmce_cold_start(GCMCoupler_ModelE *self, int yeari, int itimei, double dtsrc);

extern "C"
void gcmce_couple_native(GCMCoupler_ModelE *self,
int itime,
bool run_ice);    // if false, only initialize





}}

