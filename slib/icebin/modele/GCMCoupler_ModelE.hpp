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

struct ModelEOutputs
{
    // Pointers to arrys within ModelE

    // gcm_ovalsE[ovar](i, j, ihc)    Fortran-order 1-based indexing
    std::vector<std::unique_ptr<blitz::Array<double,3>>> gcm_ovalsE;
};

struct ModelEInputs
{
    // Pointers to arrys within ModelE

    // --------- Flux stuff
    // gcm_ivalsAI[A/E][ivar](i, j, ihc)    Fortran-order 1-based indexing
    std::vector<std::unique_ptr<blitz::Array<double,2>>> gcm_ivalsA;
    std::vector<std::unique_ptr<blitz::Array<double,3>>> gcm_ivalsE;

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
};



class DomainDecomposer_ModelE {
    ibmisc::Domain domainA_global;
    size_t ndomain;
    blitz::Array<int,1> rank_of_j;    // indexing base=1
public:

    DomainDecomposer_ModelE(std::vector<int> const &endj, ibmisc::Domain const &_domainA_global);

    /** Number of domains */
    size_t size() const { return ndomain; }

    /** Returns the MPI rank of grid cell */
    int get_domain(long ix) const {    // zero-based
        auto im_world(domainA_global[0].end);
        auto jm_world(domainA_global[1].end);

        int j = (ix / im_world) % jm_world;    // +0 for 0-based indexing
        return rank_of_j(ix);
    }
};


class GCMCoupler_ModelE : public GCMCoupler
{
public:
    double dtsrc;
    ModelEParams rdparams;    // Params straight from the rundeck (came during init)
    std::unique_ptr<boost::mpi::communicator> world;

    /** On root: separate global stuff back into individual domains.
    Works for A and E grids. */
    std::unique_ptr<DomainDecomposer_ModelE> domains;

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
    virtual ~GCMCoupler_ModelE() {}

    // Called from LISnow::allocate()
    GCMCoupler_ModelE();

    int _read_nhc_gcm();

    // The gcmce_xxx() functions do not need to be declared here
    // because everything in this class is public.


    // 1. Copies values back into modele_inputs.gcm_ivals
    void update_gcm_ivals(GCMInput const &out);

    /** Called from gcmce_Xxx() functions to update fhc, eleveE, focean, etc.
    when the ice model changes.  For now, a NOP. */
    virtual void update_topo() {}


};    // class GCMCouler_ModelE


}}

