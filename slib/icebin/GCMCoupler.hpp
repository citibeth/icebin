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

#include <cstdlib>

#include <ibmisc/VarTransformer.hpp>
#include <ibmisc/DynArray.hpp>
#include <ibmisc/udunits2.hpp>
#include <ibmisc/ConstantSet.hpp>

#include <icebin/IceCoupler.hpp>
#include <icebin/GCMParams.hpp>
#include <icebin/GCMPerIceSheetParams.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/VarSet.hpp>

namespace icebin {

struct VectorSparseParallelVectors {
    // Stores a bunch of parallel sparse vectors
    // Index of each element in the parallel vectors
    std::vector<long> index;
    // Values of for each element in the vectors.  
    std::vector<double> vals;
    // Number of _vals element per _ix element
    int nvar;
};

struct ArraySparseParallelVectorsE {
    // Stores a bunch of parallel sparse vectors
    // Index of each element in the parallel vectors
    blitz::Array<long,1> ixA;
    blitz::Array<int,1> ixHC;
    std::vector<blitz::Array<double,1>> values;
};

// extern ArraySparseParallelVectors vector_to_array(VectorSparseParallelVectors vecs);


struct GCMCouplerOutput {
    // Outputs from IceCoupler, transformed and regridded back to E/A
//    std::vector<SparseVector> gcm_ivals;    // both A and E grids.


    // Mapping from the index of a variable in gcm_ivalsE/gcm_ivalsA
    // and the index within the GCMCoupler::gcm_inputs
    VectorSparseParallelVectors[GCMCoupler::GCMI::COUNT] gcm_ivals; // gcm_ivalsE, gcm_ivalsA


    // Values required to update TOPO, etc. in ModelE
    // (see add_fhc.py for how these are to be used)
    // We can get these from AvE
    // SparseVector wAvE;     // Area of A grid cells that overlap ice
    //SparseVector areaA;    // Total (native) area of A grid cells
    //SparseVector elevA;    // Used for atmosphere orography (ZATMO in ModelE)

    // Regrid matrix to go from last step's elevation classes to this
    // step's elevation classes.
    SparseMatrix E1vE0;

    // Regrid matrix to convert to atmosphere.
    // (EvA is assumed by GCM, as long as AvE is local; see Fischer&Nowicki 2014)
    WeightedSparse AvE;

    // Used for temperature downscaling according to a lapse rate
    SparseVector elevE;
};



class GCMCoupler {
public:
    /** Type tags for subclasses of GCMCoupler */
    BOOST_ENUM_VALUES( Type, int,
        (MODELE)        (0)
        (CESM)          (1)
    );
    Type const type;

    /** Filename this coupler (including grid) was read from. */
    std::string icebin_in;

    /** Main access to the core regridding of Icebin */
    GCMRegridder regridder;

    /** Parameters (not physical constants) passed from the GCM
    through to the ice model.  These parameters cannot be specific to
    either the ice model or the GCM. */
    GCMParams gcm_params;

    /** See regridder.sheets_index */
    std::vector<std::unique_ptr<IceCoupler>> ice_couplers;

    ibmisc::UTSystem ut_system;     //!< Unit system for ConstantSets and CouplingContracts
    ibmisc::ConstantSet gcm_constants;      //!< Constants provided by the GCM

    /** Description of fields we receive from the GCM; all on the E grid. */
    VarSet gcm_outputsE;

    /** Description of fields to send back to the GCM; some on E, some on A */
    enum class GCMI { E, A, COUNT };
    VarSet[GCMI::COUNT] gcm_inputs;    // gcm_inputsE, gcm_inputsA

    /** Names of items used in the SCALARS dimension of VarTranslator.
    Used for ice_input and gcm_inputs.
    Scalars could be (eg) a timestep dt that is not known until runtime. */
    VarSet scalars;

    // Fields we read from the config file...

    GCMCoupler(Type _type) :
        type(_type),
        ut_system("")
    {
        gcm_constants.init(&ut_system);

        // Icebin requires orography on the ice grid, in order to
        // regrid in elevation space when things change.  Therefore, this
        // is added to the contract for all GCMs
        // gcm_inputs.add_field("elev2", "m", "ICE", "ice upper surface elevation");
        // No... this messes up what the GCM expects, and it's not used by the GCM.
        // Therefore, it should not be listed as a GCM input, it's a Icebin input.
        // Icebin will require it, somehow, as an IceCoupler output, and get it
        // directly from there.
    }

    virtual ~GCMCoupler() {}

    /** Read per-ice-sheet parameters that depend on the type of GCMCoupler. */
    virtual std::unique_ptr<GCMPerIceSheetParams>
    read_gcm_per_ice_sheet_params(
        ibmisc::NcIO &ncio,
        std::string const &sheet_vname) const = 0;

    virtual void ncread(
        std::string const &fname,
        std::string const &vname,
        ibmisc::Domain<int> &&domainA);

    void set_start_time(
        ibmisc::time::tm const &time_base,
        double time_start_s);

GCMCouplerOutput GCMCoupler::couple(
// Simulation time [s]
double time_s,
ArraySparseParallelVectors const &gcm_ovalsE,
bool do_run)

};

}
