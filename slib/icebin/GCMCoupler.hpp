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



struct GCMCouplerOutput {
    typedef spsparse::VectorCooArray<long, double, 2> SparseMatrix;
    typedef spsparse::VectorCooArray<long, double, 1> SparseVector;

    // http://www.boost.org/doc/libs/1_62_0/libs/serialization/doc/serialization.html#constructors
    friend class boost::serialization::access;

    // Mapping from the index of a variable in gcm_ivalsE/gcm_ivalsA
    // and the index within the GCMCoupler::gcm_inputs
    std::array<VectorMultivec<IndexAT,double>, GridAE::count> gcm_ivalsAE;

    // This is all for elevation classes in ICE space (ice_nhc, not gcm_nhc)

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
    SparseMatrix AvE1;
    SparseVector wAvE1;

    // Used for temperature downscaling according to a lapse rate
    SparseVector elevE1;

    GCMCouplerOutput(int nvar) :
        gcm_ivals({
            SparseMutlivec(nvar),
            SparseMultivec(nvar),
        }) {}
    }

    size_t nvar() { return gcm_ivals[0].nvar; }

    template<class ArchiveT>
    void serialize(ArchiveT &ar, const unsigned int file_version)
    {
        ar & gcm_ivals;
        ar & E1vE0;
        ar & AvE1;
        ar & wAvE1;
        ar & elevE1;
    }
};




struct IceOutMsgE {
    int index[3];    // i,j,ihp

    double vals[1];     // Always at least one val; but this could be extended, based on # of inputs

    double &operator[](int i) { return *(vals + i); }

    /** @return size of the struct, given a certain number of values */
    static size_t size(int nfields)
        { return sizeof(ModelEMsg) + (nfields-1) * sizeof(double); }

    static MPI_Datatype new_MPI_struct(int nfields);

    /** for use with qsort */
//  static int compar(void const * a, void const * b);

};

//struct IceCouplerOutput {
//    // Mapping between sparse and dense dimensions for E and A grids
//    std::array<SparseSet, GridAE::count> dims;
//
//    // Dense array with sparse indices represented in dims
//    // gcm_ivals[E/A](nE/nA, nvar)
//    std::array<blitz::Array<double,GridAE::count>> gcm_ivals;    // gcm_ivalsE1, gcm_ivalsA
//
//
//    // Values required to update TOPO, etc. in ModelE
//    // (see add_fhc.py for how these are to be used)
//    // We can get these from AvE
//    // SparseVector wAvE;     // Area of A grid cells that overlap ice
//    //SparseVector areaA;    // Total (native) area of A grid cells
//    //SparseVector elevA;    // Used for atmosphere orography (ZATMO in ModelE)
//
//    // Regrid matrix to go from last step's elevation classes to this
//    // step's elevation classes.
//    SparseSet dimE0;
//    EigenSparseMatrix E1vE0;
//
//    // Regrid matrix to convert to atmosphere.
//    // (EvA is assumed by GCM, as long as AvE is local; see Fischer&Nowicki 2014)
//    WeightedSparse AvE1;
//
//    // Used for temperature downscaling according to a lapse rate
//    blitz::Array<double,1> elevE1;
//};
//
//struct IceOutMsgE {
//    int index[3];    // i,j,ihp
//
//    double vals[1];     // Always at least one val; but this could be extended, based on # of inputs
//
//    double &operator[](int i) { return *(vals + i); }
//
//    /** @return size of the struct, given a certain number of values */
//    static size_t size(int nfields)
//        { return sizeof(ModelEMsg) + (nfields-1) * sizeof(double); }
//
//    static MPI_Datatype new_MPI_struct(int nfields);
//
//    /** for use with qsort */
////  static int compar(void const * a, void const * b);
//
//};




class GCMCoupler {
public:
    /** Type tags for subclasses of GCMCoupler */
    BOOST_ENUM_VALUES( Type, int,
        (MODELE)        (0)
        (CESM)          (1)
    );
    Type const type;

    // Last time this coupler was called
    double last_time_s;

    /** Filename this coupler (including grid) was read from. */
    std::string icebin_in;

    /** Main access to the core regridding of Icebin */
    GCMRegridder regridder;

    /** Parameters (not physical constants) passed from the GCM
    through to the ice model.  These parameters cannot be specific to
    either the ice model or the GCM. */
    GCMParams gcm_params;

    /** Number of elevation classes the GCM sees */
    int _nhc_gcm = -1;
    int nhc() {
        if (_nhc_gcm < 0) _nhc_gcm = get_nhc_gcm();
        return nhc_gcm();
    }
    virtual int get_nhc_gcm();

    /** See regridder.sheets_index */
    std::vector<std::unique_ptr<IceCoupler>> ice_couplers;

    ibmisc::UTSystem ut_system;     //!< Unit system for ConstantSets and CouplingContracts
    ibmisc::ConstantSet gcm_constants;      //!< Constants provided by the GCM

    /** Description of fields we receive from the GCM; all on the E grid. */
    VarSet gcm_outputsE;

    /** Description of fields to send back to the GCM; some on E, some on A */
    std::array<VarSet, GridAE::count> gcm_inputs;    // gcm_inputsE, gcm_inputsA

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

void GCMCoupler::couple(
// Simulation time [s]
double time_s,
ArraySparseParallelVectorsE const &gcm_ovalsE,
MultiDomainGCMCoupleOutput &out,
bool do_run)

};

std::vector<GCMCouplerOutput> split_by_domain(
    GCMCouplerOutput const &out,
    DomainDecomposer const &domainsA,
    DomainDecomposer const &domainsE);

}
