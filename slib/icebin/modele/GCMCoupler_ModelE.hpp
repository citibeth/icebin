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

#include <icebin/GCMCoupler.hpp>

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


class GCMPerIceSheetParams_ModelE : public icebin::GCMPerIceSheetParams {
public:
    ModelE_CouplingType coupling_type;
};
#endif


// Parameters read from the ModelE rundeck
static int const MAX_CHAR_LEN = 128;    // From Dictionary_mod.F90
struct ModelEParams
{
    char icebin_segments[MAX_CHAR_LEN];
    double dtsrc;
    int yeari;
};

struct ModelEOutputs
{
    // Pointers to arrys within ModelE

    // gcm_ovalsE[ovar](i, j, ihc)    Fortran-order 1-based indexing
    std::vector<std::unique_ptr<blitz::Array<double,3>>> gcm_ovalsE;
};

struct ModelEInputs
{
    std::vector<HCSegmentData> hc_segments;
    SegmentData *icebin_segment;    // The segment with IceBin-generated elevation classes

    // Relationship between elevation classes in ModelE and elevation classes in IceBin
    int const icebin_base_hc;    // First GCM elevation class that is an IceBin class (0-based indexing)
    int const nhc_gcm;    // Number of elevation classes used by ModelE

    // Pointers to arrys within ModelE

    // --------- Flux stuff
    // gcm_ivalsAI[A/E][ivar](i, j, ihc)    Fortran-order 1-based indexing
    std::vector<std::vector<std::unique_ptr<blitz::Array<double,3>>>> gcm_ivals;

    // --------- State variables we can update inside ModelE
    // i,j,ihc arrays on Elevation grid
    blitz::Array<double,3> fhc;
    blitz::Array<double,3> elevE;

    // i,j arrays on Atmosphere grid
    blitz::Array<double,2> focean;
    blitz::Array<double,2> flake;
    blitz::Array<double,2> fgrnd;    // Alt: fearth0
    blitz::Array<double,2> fgice;    // Alt: flice
    blitz::Array<double,2> zatmo;      // i,j

    void ModelEInputs::update_gcm_ivals(GCMCouplerOutput const &out);
};



class DomainDecomposer_ModelE {
    size_t ndomain;
    blitz::Array<int,1> rank_of_j;    // indexing base=1
public:

    DomainDecomposer_ModelE(std::vector<int> const &startj, im_world, jm_world) :    // Starts from ModelE; j indexing base=1
        rank_of_j(startj[startj.size()-1], blitz::fortranArray),    // 1-based indexing
        im_world(_im_world), jm_world(_jm_world);

    /** Number of domains */
    size_t size() { return ndomain; }

    /** Returns the MPI rank of grid cell */
    int get_rank(long ix) {    // zero-based
        int j = (ix / im_world) % jm_world;    // +0 for 0-based indexing
        return rank_of_j(ix);
    }
}


class GCMCoupler_ModelE : public GCMCoupler
{
public:
    std::unique_ptr<boost::mpi::communicator> world;

    ModelEOutputs modele_outputs;

    // Variables borrowed from ModelE, used to return data to it.
    // All these variables are Fortran-order, 1-based indexing
    ModelEInputs modele_inputs;

    // Low and high indices for this MPI rank.
    // Indices are in Fortran order (im, jm) with zero-based indexing
    ibmisc::Domain domainA;
    // Low and high indices for global domain (Fortran order, 0-based)
    ibmisc::Domain domainA_global;


public:
    // Called from LISnow::allocate()
    GCMCoupler_ModelE();

    // The gcmce_xxx() functions do not need to be declared here
    // because everything in this class is public.

#if 0
    /** Read per-ice-sheet parameters that depend on the type of GCMCoupler. */
    std::unique_ptr<GCMPerIceSheetParams>
    read_gcm_per_ice_sheet_params(
        ibmisc::NcIO &ncio,
        std::string const &sheet_vname);
#endif
};    // class GCMCouler_ModelE


}}

