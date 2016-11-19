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
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <icebin/GCMCoupler.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/contracts/contracts.hpp>
#include <spsparse/multiply_sparse.hpp>
#include <spsparse/sort.hpp>

#ifdef USE_PISM
#include <icebin/pism/IceModel_PISM.hpp>
#endif

using namespace spsparse;
using namespace ibmisc;
using namespace netCDF;

namespace icebin {

std::unique_ptr<IceCoupler> new_ice_coupler(IceCoupler::Type type,
    GCMCoupler const *_coupler, IceRegridder *_sheet)
{
    std::unique_ptr<IceCoupler> ice_coupler;

    switch(type.index()) {
#if 0
        case IceCoupler::Type::DISMAL :
            ice_coupler.reset(new IceCoupler_DISMAL);
        break;
#endif
        case IceCoupler::Type::WRITER :
            ice_coupler.reset(new IceCoupler_Writer);
        break;
#ifdef USE_PISM
        case IceCoupler::Type::PISM :
            ice_coupler.reset(new gpism::IceCoupler_PISM);
        break;
#endif
        default :
            (*icebin_error)(-1,
                "Unknown IceCoupler::Type %s", type.str());
    }


    // Do basic initialization...
    ice_coupler->coupler = _coupler;
    ice_coupler->sheet = _sheet;
    ice_coupler->ice_constants.init(&_coupler->ut_system);

    return ice_coupler;


}

std::unique_ptr<IceCoupler> new_ice_coupler(NcIO &ncio, std::string vname,
    GCMCoupler const *_coupler, IceRegridder *_sheet)
{
    std::string vn(vname + ".info");
    auto info_v = get_or_add_var(ncio, vn, "int64", {});

    IceCoupler::Type type;
    get_or_put_att_enum(info_v, ncio.rw, "ice_coupler", type);

    return new_ice_coupler(type, _coupler, _sheet);
}

void IceCoupler::ncread(
ibmisc::NcIO &ncio, std::string const &vname_sheet)
{
    gcm_per_ice_sheet_params = coupler->read_gcm_per_ice_sheet_params(ncio, vname_sheet);
}



IceCoupler::~IceCoupler() {}

#if 0
VarSet *IceCoupler::new_VarSet() {
    _extra_contracts.push_back(
        std::unique_ptr<VarSet>(
        new VarSet()));
    return _extra_contracts[_extra_contracts.size()-1].get();
}
#endif

// ==========================================================
bool IceCoupler::am_i_root() const
    { return coupler->am_i_root(); }

/** Allocate vectors in preparation of calling an ice model. */
void IceCoupler::allocate_ice_ovals_I()
{
    // Check for program errors
    if (!coupler->am_i_root()) (*icebin_error)(-1,
        "IceCoupler::allocate_ice_ovals_I() should only be called from GCM root MPI node.  Fix the code.");

    if (ice_ovals_I.size() != 0) (*icebin_error)(-1,
        "[%d] IceCoupler::allocate_ice_ovals_I(): called twice without a free() inbetween.  Fix the code. (old size is %ld)\n", coupler->gcm_params.gcm_rank, ice_ovals_I.size());

    // Allocate for direct output from ice model
    VarSet const &ocontract(contract[IceCoupler::OUTPUT]);
    for (int i=0; i < ocontract.size(); ++i)
        ice_ovals_I.push_back(blitz::Array<double,1>(ndata()));
}


/** Allocate in preparation of var transformations (but not regridding yet) */
void IceCoupler::allocate_gcm_ivals_I()
{
    // Check for program errors
    if (!coupler->am_i_root()) (*icebin_error)(-1,
        "IceCoupler::allocate_ice_ivals_I() should only be called from GCM root MPI node.  Fix the code.\n");

    if (gcm_ivals_I.size() != 0) (*icebin_error)(-1,
        "IceCoupler::allocate_gcm_ivals_I(): called twice without a free() inbetween.  Fix the code.\n");


    VarSet const &gcm_inputs(coupler->gcm_inputs);
    for (int i=0; i < gcm_inputs.size(); ++i) {
        gcm_ivals_I.push_back(blitz::Array<double,1>(ndata()));
    }
}

/** Free portions not needed after finished calling ice model and
applying variable transform.  This will be variables desired on
anything other than the ELEVATION grid. */
void IceCoupler::free_ice_ovals_I()
{
    // Check for program errors
    if (!coupler->am_i_root()) (*icebin_error)(-1,
        "IceCoupler::free_ice_ovals_I() should only be called from GCM root MPI node.  Fix the code.\n");

    ice_ovals_I.clear();
}

/** Free all memory used by this.  Called when we're done with a coupling timestep. */
void IceCoupler::free_ovals_ivals_I()
{
    // Check for program errors
    if (!coupler->am_i_root()) (*icebin_error)(-1,
        "IceCoupler::free_ovals_ovals_I() should only be called from GCM root MPI node.  Fix the code.\n");

    ice_ovals_I.clear();
    gcm_ivals_I.clear();
}

// -----------------------------------------------------------
/** Allocates and sets gcm_ivals_I variable
@param mask Control which GCM input variables to set (according to, etc, INITIAL flag in contract) */
void IceCoupler::set_gcm_inputs(unsigned int mask)
{
  printf("BEGIN IceCoupler::set_gcm_inputs()\n");
    allocate_gcm_ivals_I();

    // Compute the variable transformation
    ibmisc::VarTransformer &vt(var_transformer[IceCoupler::OUTPUT]);
    ibmisc::CSRAndUnits trans = vt.apply_scalars({
        std::make_pair("unit", 1.0)});

    VarSet const &gcm_inputs(coupler->gcm_inputs);

    // Apply the variable transformation
    for (int xi=0; xi<vt.dim(ibmisc::VarTransformer::OUTPUTS).size(); ++xi) {   // xi is index of output variable
        gcm_ivals_I[xi] = 0;    // Vector operation: clear before sum
        VarMeta const &cf(gcm_inputs[xi]);

        if ((cf.flags & mask) != mask) continue;

        // Consider each output variable separately...
        std::vector<std::pair<int, double>> const &row(trans.mat[xi]);
        for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
            int xj = xjj->first;        // Index of input variable
            double io_val = xjj->second;    // Amount to multiply it by
            gcm_ivals_I[xi] += ice_ovals_I[xj] * io_val;        // blitz++ vector operation
        }
        gcm_ivals_I[xi] += trans.units[xi];
    }

    printf("END IceCoupler::set_gcm_inputs()\n");
}
// ==========================================================
static double const nan = std::numeric_limits<double>::quiet_NaN();

// REMEMBER: Decoding converts a set of (index, value) pairs into
// normal arrays (with NaN where no value was given.)
void IceCoupler_Decode::run_timestep(double time_s,
	blitz::Array<int,1> const &indices,
	std::vector<blitz::Array<double,1>> const &ivals2)
{
printf("BEGIN IceCoupler_Decode::run_timestep(time_s = %f) size=%ld\n", time_s, indices.size());

	blitz::Array<int,1> nindices;
	std::vector<blitz::Array<double,1>> nivals2;

	blitz::Array<int,1> const *xindices;
	std::vector<blitz::Array<double,1>> const *xivals2;

	// Test out freestanding sorting/consolidation code
	// This section of code is NOT necessary.
	if (false) {
		std::vector<int> perm = sorted_perm(indices);
		nindices.reference(blitz::Array<int,1>(indices.size()));
		int nconsolidated = consolidate_by_perm(indices, perm, indices, nindices, DuplicatePolicy::REPLACE);
printf("IceCoupler_Decode: consolidated from %d down to %d\n", indices.size(), nconsolidated);
		nindices.reference(blitz::Array<int,1>(nconsolidated));
		consolidate_by_perm(indices, perm, indices, nindices, DuplicatePolicy::REPLACE);

		for (unsigned int i=0; i<ivals2.size(); ++i) {
			blitz::Array<double,1> nvals(nconsolidated);
			consolidate_by_perm(indices, perm, ivals2[i], nvals, DuplicatePolicy::ADD);
			nivals2.push_back(nvals);
		}

		xindices = &nindices;
		xivals2 = &nivals2;
	} else {
		xindices = &indices;
		xivals2 = &ivals2;
	}

	std::vector<blitz::Array<double,1>> ivals2d;	/// Decoded fields

	// Naming convention on array variables:
	//     ivals2 = Vector of Values-arrays on grid2 (ice grid)
	//     ivals2d = Vector of DECODED values-arrays on grid2
	//     vals = Individual value array from ivals2
	//     valsd = Individual valu array from ivals2d
	// Loop through the fields we require
	VarSet const &icontract(contract[IceCoupler::INPUT]);
	for (int i=0; i<icontract.size(); ++i) {

		blitz::Array<double,1> const &vals((*xivals2)[i]);

		// Decode the field!
		blitz::Array<double,1> valsd(ndata());
		valsd = nan;
		int n = xindices->size();
		for (int i=0; i < n; ++i) {
			int ix = (*xindices)(i);
			// Do our own bounds checking!
			if (ix < 0 || ix >= ndata()) (*icebin_error)(-1,
                "IceCoupler: index %d out of range [0, %d)\n", ix, ndata());

			// Add this value to existing field
			double &oval = valsd(ix);
			if (std::isnan(oval)) oval = vals(i);
			else oval += vals(i);
		}

		// Convert any remaining nans to default value,
		// so we have a valid number everywhere.
		double default_value = icontract[i].default_value;
		for (int j=0; j<ndata(); ++j) {
			double &val(valsd(j));
			if (std::isnan(val)) val = default_value;
		}

		// Store decoded field in our output
		ivals2d.push_back(valsd);
printf("Done decoding required field, %s\n", icontract[i].name.c_str());
	}

	// Pass decoded fields on to subclass
	run_decoded(time_s, ivals2d);
printf("END IceCoupler_Decode::run_timestep(%f)\n", time_s);
}

// ------------------------------------------------------------

SparseMatrix regrid_M(IceRegridder &regridder, std::string const &spec_name)
    { return std::move(regridder.regrid(spec_name)->M); }



// ==============================================================
/** @return ice_ovalsI */
void IceCoupler::couple(
double time_s,
// Values from GCM, passed GCM -> Ice
blitz::Array<long,1> gcm_ovalsE_index,    // Indices of sparse vectors...
std::vector<blitz::Array<double,1>> gcm_ovalsE_values,    // values[var](i)
GCMCoupleOutput &out)    // Accumulate matrices here...
{
    // Store regridding matrices for the last timestep, which we will
    // need to create ice_ivals
//    std::unique_ptr<WeightedSparse> IvEd0(std::move(IvEd));    // Dense A and E
//    std::unique_ptr<WeightedSparse> IvAd0(std::move(IvAd));

    // ========== Get Ice Inputs
    blitz::Array<double,2> ice_ivalsI({ndata(), contract[INPUT].size()});
    blitz::Array<double,2> ice_ovalsI({ndata(), contract[OUTPUT].size()});
    ice_ivalsI = 0;
    ice_ovalsI = 0;

    // Save old IvE for now
    auto IvE0(std::move(IvE));

    // ------------- Form ice_ivalsI
    // Densify gcm_ovals and IvE
    SparseSetT dimE;
    SparseMatrix IvEd0;
    dimE.add_sorted(IvE0.dim_begin(1), IvE0.dim_end(1));
    densify_one_dim(IvEd0, IvE0, dimE, 1);

    // Densify gcm_ovalsE --> gcm_ovalsEd
    blitz::Array<double,2> gcm_ovalsEd(dimE.sparse_extent(), gcm_cupler->gcm_outputsE);
    gcm_ovalsEd = 0;
    for (size_t i=0; i<gcm_ovalsE_index.size(); ++i) {
        iEd = dimE.to_dense(gcm_ovalsE_index[i]);
        for (size_t ivar=0; ivar<gcm_ovalsE_values.size(); ++i) {
            gcm_ovals(iEd, ivar) += gcm_ovalsE_values[ivar](i);
        }
    }

    // Get the CSR sparse matrix to convert GCM outputs to ice model inputs
    auto scalars({
        std::make_pair("by_dt", 1.0 / ((itime - api->itime_last) * api->dtsrc)),
        std::make_pair("unit", 1.0)});
    CSRAndUnits icei_v_gcmo(var_transformer[INPUT].apply_scalars(scalars));

    // Regrid & combine to form ice_ivalsI
    for (auto iM = IvEd.begin(); iM != IvEd.end(); ++iM) {
        // Transform units on the input while multiplying by M
        for (size_t ovar = 0; ovar < contracts[INPUT].size(); ++ovar) {
            double zval = 0;
            std::vector<std::pair<int, double>> const &row(icei_v_gcmo.mat[xi]);
            for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
                int xj = xjj->first;
                double xval = xjj->second;
                zval += xval * gcm_ovalsEd(iM.index(1), xj);
            }
            ice_ivalsI(iM.index(0), ovar) += iM.val() * zval;
        }
    }

    // ========= Step the ice model forward
    iwriter()->run_timestep(time_s, ice_ivalsI);    // TODO: Change to a write that takes a contract and a vector of fields
    run_timestep(time_s, ice_ivalsI, ice_ovalsI);
    owriter()->run_timestep(time_s, ice_ovalsI);

    // ========== Update regridding matrices
    _update_elevI();
    RegridMatrices rm(regridder);

    // Compute IvE (for next timestep)
    IvE = regrid_M(re, "IvE");
    dimE.add_sorted(IvE->M.dim_begin(1), IvE->M.dim_end(1));
    densify_one_dim(IvEd, IvE->M, dimE, 1);

    // Compute regrid matrices we need now
    auto EvI(regrid_M(re, "EvI"));
    auto AvI(regrid_M(re, "AvI"));
    auto AvE(regrid_M(re, "AvE"));

    // Accumulate global things for all ice sheets....
    multiply(coupler_ret.E1vE0, EvI, IvE0);
    copy(coupler_ret.AvE.M, AvE.M);
    copy(coupler_ret.AvE.weight, AvE.weight);
    multiply(coupler_ret.elevE, EvI, elevI);


    // ========= Compute gcm_ivalsE = EvI * vt * ice_ovals
    CSRAndUnits gcmi_v_iceo(var_transformer[OUTPUT].apply_scalars(scalars));

    std::array<SparseMatrix *, GCMCoupler::GCMI::COUNT> XvIs;
    XvIs[GCMCoupler::GCMI::E] = &EvI;
    XvIs[GCMCoupler::GCMI::A] = &AvI;
    // Do it once for _E variables and once for _A variables.
    for (int i=0; i < GCMCoupler::GCMI::COUNT; ++i) {
        VarSet &contract(gcm_coupler->gcm_inputs[i]);
        SparseMatrix &XvI(*XvIs[i]);
        SparseParallelVectors &gcm_ivalsX(coupler_ret.gcm_ivals[i]);

        // Do the multiplication
        for (auto iM = XvI.begin(); iM != XvI.end(); ++iM) {

            // Transform units on the input while multiplying by M
            gcm_ivalsX.index.push_back(iM.index(0));
            for (size_t ovar = 0; ovar < contract.size(); ++ovar) {
                double zval = 0;
                std::vector<std::pair<int, double>> const &row(gcmi_v_iceo.mat[xi]);
                for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
                    int xj = xjj->first;
                    double xval = xjj->second;
                    zval += xval * ice_ovalsI(iM.index(1), xj);
                }
                gcm_ivalsX.vals.push_back(iM.val() * zval);
            }
        }
    }
}
