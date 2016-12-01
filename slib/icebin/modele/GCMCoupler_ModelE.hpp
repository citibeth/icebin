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

struct ModelEOutputs
{
    // Pointers to arrys within ModelE

    // gcm_ovalsE[ovar](i, j, ihc)    Fortran-order 1-based indexing
    std::vector<std::unique_ptr<blitz::Array<double,3>>> gcm_ovalsE;
};

struct ModelEInputs
{
    // Pointers to arrys within ModelE

    // gcm_ivalsAI[A/E][ivar](i, j, ihc)    Fortran-order 1-based indexing
    std::vector<std::vector<std::unique_ptr<blitz::Array<double,3>>>> gcm_ivals;

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




class GCMCoupler_ModelE : public GCMCoupler
{
    ModelEOutputs modele_outputs;

    // Variables borrowed from ModelE, used to return data to it.
    // All these variables are Fortran-order, 1-based indexing
    ModelEInputs modele_inputs;

    // The first GCM elevation class that is an IceBin class (0-based indexing)
    int icebin_base_hc;
    int icebin_nhc;    // Number of elevation classes used by IceBin

    // Low and high indices for this MPI node.
    // Indices are in Fortran order (im, jm) with zero-based indexing
    ibmisc::Domain<int> domainA;
    // Low and high indices for global domain (Fortran order, 0-based)
    ibmisc::Domain<int> domainA_global;

public:
    GCMCoupler_ModelE();

    /** Read per-ice-sheet parameters that depend on the type of GCMCoupler. */
    std::unique_ptr<GCMPerIceSheetParams>
    read_gcm_per_ice_sheet_params(
        ibmisc::NcIO &ncio,
        std::string const &sheet_vname);
};


}}

