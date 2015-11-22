/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <array>
#include <boost/enum.hpp>
#include <boost/filesystem/path.hpp>
#include <giss/time.hpp>
#include <giss/VarTransformer.hpp>
#include <glint2/IceSheet.hpp>
#include <giss/CouplingContract.hpp>
#include <giss/ConstantSet.hpp>
#include <glint2/GCMPerIceSheetParams.hpp>

namespace glint2 {

class GCMCoupler;
class IceModel_Writer;

// -------------------------------------------------------
// GCM-specific types for parameters to setup_contracts_xxxxx()

namespace modele{
	class GCMCoupler_ModelE;
	class GCMPerIceSheetParams_ModelE;
}
// -------------------------------------------------------
class IceModel {

friend std::unique_ptr<IceModel> read_icemodel(
	std::string const &name,
	GCMCoupler const *coupler,
	NcFile &nc,
	std::string const &vname,
	std::unique_ptr<GCMPerIceSheetParams> &&gcm_per_ice_sheet_params,
	IceSheet *sheet);

public:
	/** The name given to this IceModel, used to index */
	std::string const name;

protected:

	/** The grid for this IceModel. */
	std::shared_ptr<glint2::Grid> grid2;

	// Parameters provided by the GCM, to inform the coupling
	std::unique_ptr<GCMPerIceSheetParams> gcm_per_ice_sheet_params;
	GCMCoupler const * const coupler;		// parent back-pointer

	// Writers called to record the input and output seen by this IceModel
	std::unique_ptr<IceModel_Writer> _iwriter, _owriter;

public:
	BOOST_ENUM_VALUES( Type, int,
		(DISMAL)		(0)		// Demo Ice Sheet Model and LandIce
		(PISM)			(1)
		(ISSM)			(2)
		(WRITER)		(3)
	);
	const IceModel::Type type;

	enum IO {INPUT, OUTPUT};

	/** Constants obtained from the GCM */
	giss::ConstantSet ice_constants;

	/** Ordered specification of the variables (w/ units)
	to be passed GLINT2->IceModel */
	std::array<giss::CouplingContract, 2> contract;		// [INPUT|OUTPUT]
	std::array<giss::VarTransformer, 2> var_transformer;

	/** Placeholder for additional coupling contracts that had to be allocated. */
	std::vector<std::unique_ptr<giss::CouplingContract>> _extra_contracts;

	// --------------------------------------------
	// Buffers used to receive ice model output, and regrid it.

	/** Direct output from the ice model (on the ice grid).
	This is allocated (for the ice model ROOT node only) inside
	run_timestep() */
	std::vector<blitz::Array<double,1>> ice_ovals_I;

	/** Input to the GCM, but on the ice grid */
	std::vector<blitz::Array<double,1>> gcm_ivals_I;
	// ======================================================

	void set_writers(
		std::unique_ptr<IceModel_Writer> &&iwriter,
		std::unique_ptr<IceModel_Writer> &&owriter);

	/** @return the iwriter associated with this IceModel, if the
	IceModel is NOT responsible for calling the writer itself.
	If the IceModel WILL call the writer, returns NULL. */
	virtual IceModel_Writer *iwriter() { return _iwriter.get(); }
	virtual IceModel_Writer *owriter() { return _owriter.get(); }


	/** Tells whether we are running on the root node of the ice model. */
	virtual bool am_i_root() const;

	/** Allocate vectors in preparation of calling an ice model (ROOT only). */
	void allocate_ice_ovals_I();

	/** Allocate in preparation of var transformations (but not regridding yet) (ROOT only) */
	void allocate_gcm_ivals_I();

	/** Free portions not needed after finished calling ice model and
	applying variable transform.  This will be variables desired on
	anything other than the ELEVATION grid. (ROOT only) */
	void free_ice_ovals_I();

	/** Free all memory used by this.  Called when we're done with a coupling timestep. (ROOT only) */
	void free_ovals_ivals_I();

	/** Allocates and sets gcm_ivals_I variable */
	void set_gcm_inputs(unsigned int mask);

	// --------------------------------------------
	/** Allocate a new giss::CouplingContract, with the same lifetime as this IceModel. */
	giss::CouplingContract *new_CouplingContract();

	IceModel(IceModel::Type _type, std::string const &_name, GCMCoupler const *_coupler);
	virtual ~IceModel();

	long ndata() const { return grid2->ndata(); }

	// --------------------------------------------------
	// GCM-specific methods used to set up the contract for
	// a particular GCM-IceModel pair
	virtual void setup_contracts_modele();

	// --------------------------------------------------

protected:
	/** Call this from main init() method (below) */
	void init(std::shared_ptr<glint2::Grid> const &_grid2)
		{ this->grid2 = _grid2; }
public:

	/** Initialize any grid information, etc. from the IceSheet struct.
	@param vname_base Construct variable name from this, out of which to pull parameters from netCDF */
	virtual void init(
//		IceModel::GCMParams const &gcm_params,
		std::shared_ptr<glint2::Grid> const &grid2,
		NcFile &nc,
		std::string const &vname_base)
	{
		fprintf(stderr, "IceModel::init() must be implemented!\n");
		giss::exit(1);
	}

	/** Event handler to let IceModels know the start time is (finally) set */
	virtual void start_time_set() {}

	/** Run the ice model for one coupling timestep.
	@param time_s Seconds since GCMParams::time_base.  Helps with debugging.
	@param index Index of each input grid value in ivals2.
	@param ivals2 The values themselves (sparse representation).
           Their meaning (SMB, T, etc) is determined
           by the place in the array, as specified by the appropriate
           INPUT contract for this ice model.
	*/
	virtual void run_timestep(double time_s,
		blitz::Array<int,1> const &indices,
		std::vector<blitz::Array<double,1>> const &ivals2) = 0;

	/** Called at the beginning.  Returns variables in same place as run_timestep(),
	but doesn't actually run the timestep. */
	virtual void get_initial_state(double time_s) {
		fprintf(stderr, "get_initial_state() not implemented.\n");
		giss::exit(1);
	}


	/** Allows the IceModel to change the inputs used to create the
	regridding transformations.  This is used, for example, to make
	elev2 and mask2 consistent with an existing ice model (eg, PISM).
	It is called after init().
	Default implementation is to do nothing. */
	virtual void update_ice_sheet(NcFile &nc, std::string const &vname,
		IceSheet *sheet) {}

};

extern std::unique_ptr<IceModel> read_icemodel(
	std::string const &name,
	GCMCoupler const *coupler,
	NcFile &nc, std::string const &vname,
	std::unique_ptr<GCMPerIceSheetParams> &&gcm_per_ice_sheet_params,
	IceSheet *sheet = NULL);

}
