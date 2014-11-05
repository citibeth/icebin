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

#include <cstdlib>
#include <giss/DynArray.hpp>
#include <giss/Dict.hpp>
#include <giss/udunits2.hpp>
#include <glint2/IceModel.hpp>
#include <glint2/IceModel_Writer.hpp>
#include <giss/ConstantSet.hpp>
#include <glint2/GCMParams.hpp>
#include <glint2/GCMPerIceSheetParams.hpp>
#include <glint2/GridDomain.hpp>
#include <glint2/MatrixMaker.hpp>

namespace glint2 {




struct SMBMsg {
	int sheetno;
	int i2;			// Index into ice model
	double vals[1];		// Always at least one val; but this could be extended

	double &operator[](int i) { return *(vals + i); }

	/** @return size of the struct, given a certain number of values */
	static size_t size(int nfields)
		{ return sizeof(SMBMsg) + (nfields-1) * sizeof(double); }

	static MPI_Datatype new_MPI_struct(int nfields);

	/** for use with qsort */
	static int compar(void const * a, void const * b);

};


class GCMCoupler {
public:
	/** Type tags for subclasses of GCMCoupler */
	BOOST_ENUM_VALUES( Type, int,
		(MODELE)		(0)
		(CESM)			(1)
	);
	Type const type;

	/** Main access to the core regridding of Glint2 */
	std::unique_ptr<MatrixMaker> maker;

	/** Parameters passed from the GCM through to the ice model.
	These parameters cannot be specific to either the ice model or the GCM. */
	GCMParams gcm_params;
	giss::MapDict<int,IceModel> models;

	/** Associated data structures to write out the exact inputs seen
	by each ice model. */
	giss::MapDict<int,IceModel_Writer> writers[2];	// INPUT and OUTPUT of ice model

	giss::UTSystem ut_system;		//!< Unit system for ConstantSets and CouplingContracts
	giss::ConstantSet gcm_constants;		//!< Constants provided by the GCM

	/** Fields we receive from the GCM */
	giss::CouplingContract gcm_outputs;

	/** Fields to send back to the GCM */
	giss::CouplingContract gcm_inputs;

	/** Names of items used in the SCALARS dimension of VarTranslator.
	Used to convert GCM outputs in terms of (eg) coupling timestep (known to the GCM) */
	giss::CouplingContract ice_input_scalars;

// Not needed, we use ice_input_scalars instead.
//	/** Names of items used in the SCALARS dimension of VarTranslator.
//	Used to convert ice outputs to GCM inputs in terms of (eg)
//	coupling timestep (known to the GCM) */
//	giss::CouplingContract gcm_input_scalars;

	// Fields we read from the config file...

	/** File to which to write gcm_output.  If "", then don't write. */
	std::string gcm_out_file;

	GCMCoupler(Type _type) :
		type(_type), ut_system(""),
		gcm_constants(&ut_system)
	{

		// Glint2 requires orography on the ice grid, in order to
		// regrid in elevation space when things change.  Therefore, this
		// is added to the contract for all GCMs
		gcm_inputs.add_field("elev2", "m", "ICE", "ice upper surface elevation");
	}

	virtual ~GCMCoupler() {}

	/** Read per-ice-sheet parameters that depend on the type of GCMCoupler. */
	virtual std::unique_ptr<GCMPerIceSheetParams>
	read_gcm_per_ice_sheet_params(
		NcFile &nc,
		std::string const &sheet_vname) = 0;

	/** Called from just after read_from_netcdf().  GCM-specific implementation
	sets up the contracts for this GCM - IceModel pair.  (Usually by calling
	through to GCM-specific virtual methods on the IceModel. */
	virtual void setup_contracts(IceModel &ice_model) const = 0;

	virtual void read_from_netcdf(
		NcFile &nc, std::string const &vname,
		std::unique_ptr<GridDomain> &&mdomain);

	void set_start_time(
		giss::time::tm const &time_base,
		double time_start_s);


	/** Returns a unique rank number for each node in the parallel computation.
	Useful for debugging-type output. */
	int rank() const;

protected:
	/** @param time_s Time since start of simulation, in seconds
	Fields contained in the SMBMsg are INPUTS from the GCM.  They are
	therefore arranged according to gmc_inputs.  GCM inputs are converted
	into ice model inputs within IceModel::run_timestep(), which
	is called at the end of this method.
	@see gmc_inputs*/
	void call_ice_model(
		IceModel *model,
		int sheetno,
		double time_s,
		giss::DynArray<SMBMsg> &rbuf,
		SMBMsg *begin, SMBMsg *end);


public:
	/** Top-level general-purpose method, called by glint2_modele.cpp
	(or other GCM-specific code).
	@param time_s Time (seconds) since the start of the GCM run.
	@param nfields Number of fields in sbuf.  Not all will necessarily be filled, in the case of heterogeneous ice models.
	@param sbuf the (filled) array of ice grid values for this MPI node.
	*/
	void couple_to_ice(double time_s,
		int nfields,
		giss::DynArray<SMBMsg> &sbuf,
		std::vector<blitz::Array<double,1>> &gcm_ivals);
};

}
