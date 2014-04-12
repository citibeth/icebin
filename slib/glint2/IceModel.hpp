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

#include <boost/enum.hpp>
#include <boost/filesystem/path.hpp>
#include <giss/time.hpp>
#include <giss/VarTransformer.hpp>
#include <glint2/IceSheet.hpp>
#include <glint2/CouplingContract.hpp>

namespace glint2 {

// -------------------------------------------------------
// GCM-specific types for parameters to setup_contract_xxxxx()

namespace modele{
	class GCMCoupler_ModelE;
	class ContractParams_ModelE;
}
// -------------------------------------------------------
class IceModel {
public:
	/** Parameters passed from the GCM through to the ice model.
	These parameters cannot be specific to either the ice model or the GCM. */
	struct GCMParams {
		MPI_Comm gcm_comm;		// MPI communicator used by the GCM
		int gcm_rank;			// Rank of this process in gcm_comm
		int gcm_root;			// Rank of root process in gcm_comm
		boost::filesystem::path config_dir;	// Where to look for Ice Model configuration files
		giss::time::tm time_base;	// Corresponds to time_s == 0
		double time_start_s;		// Start of simulation, as far as ice model is concerned (seconds since time_base).

		GCMParams();

		GCMParams(
			MPI_Comm const _gcm_comm,
			int _gcm_root,
			boost::filesystem::path const &_config_dir,
			giss::time::tm const &_time_base,
			double time_start_s);
	};

protected:
	GCMParams gcm_params;
	
public:
	BOOST_ENUM_VALUES( Type, int,
		(DISMAL)		(0)		// Demo Ice Sheet Model and LandIce
		(PISM)			(1)
		(ISSM)			(2)
	);
	const IceModel::Type type;


	enum IO {INPUT, OUTPUT};

	/** Ordered specification of the variables (w/ units)
	to be passed GLINT2->IceModel */
	CouplingContract contract[2];		// [INPUT|OUTPUT]
	giss::VarTransformer var_transformer[2];

	/** Placeholder for additional coupling contracts that had to be allocated. */
	std::vector<std::unique_ptr<CouplingContract>> _extra_contracts;

	/** Allocate a new CouplingContract, with the same lifetime as this IceModel. */
	CouplingContract *new_CouplingContract() {
		_extra_contracts.push_back(
			std::unique_ptr<CouplingContract>(
			new CouplingContract));
		return _extra_contracts[_extra_contracts.size()-1].get();
	}

	IceModel(IceModel::Type _type) : type(_type) {}
	virtual ~IceModel() {}

	// --------------------------------------------------
	// GCM-specific methods used to set up the contract for
	// a particular GCM-IceModel pair
	virtual void setup_contract_modele(
		glint2::modele::GCMCoupler_ModelE const &coupler,
		glint2::modele::ContractParams_ModelE const &params)
	{
		fprintf(stderr, "setup_contract_modele() not implemented for IceModel type %s", type.str());
		throw std::exception();
	}


	/** Non-GCM-specific portion of setting up contracts.  Called by
	GCMCoupler_yyy setup AFTER setup_contract_yyy(). */
	virtual void finish_contract_setup()
		{ printf("IceModel::finish_contract_setup() does nothing.\n"); }

	// --------------------------------------------------

	void init(IceModel::GCMParams const &_gcm_params)
	{ gcm_params = _gcm_params; }

	/** Initialize any grid information, etc. from the IceSheet struct.
	@param vname_base Construct variable name from this, out of which to pull parameters from netCDF */
	virtual void init(
		IceModel::GCMParams const &gcm_params,
		std::shared_ptr<glint2::Grid> const &grid2,
		NcFile &nc,
		std::string const &vname_base,
		NcVar *const_var)
	{
		fprintf(stderr, "IceModel::init() must be implemented!\n");
		throw std::exception();
	}

	/** @param index Index of each grid value.
	@param vals The values themselves -- could be SMB, Energy, something else...
	TODO: More params need to be added.  Time, return values, etc.
	@param time_s Seconds since GCMParams::time_base
	Helps with debugging. */
	virtual void run_timestep(double time_s,
		blitz::Array<int,1> const &indices,
		std::vector<blitz::Array<double,1>> const &vals2) = 0;

	/** Allows the IceModel to change the inputs used to create the
	regridding transformations.  This is used, for example, to make
	elev2 and mask2 consistent with an existing ice model (eg, PISM).
	It is called after init().
	Default implementation is to do nothing. */
	virtual void update_ice_sheet(NcFile &nc, std::string const &vname,
		IceSheet *sheet) {}

};

extern std::unique_ptr<IceModel> read_icemodel(
	IceModel::GCMParams const &gcm_params,
	NcFile &nc, std::string const &vname, IceSheet *sheet = NULL);

}
