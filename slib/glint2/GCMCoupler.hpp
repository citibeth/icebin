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
#include <boost/filesystem.hpp>
#include <giss/DynArray.hpp>
#include <giss/Dict.hpp>
#include <glint2/IceModel.hpp>
#include <glint2/CouplingContract.hpp>

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

	IceModel::GCMParams const gcm_params;
	giss::MapDict<int,IceModel> models;

	/** Fields we receive from the GCM */
	CouplingContract gcm_outputs;

	GCMCoupler(Type _type, IceModel::GCMParams const &_gcm_params) :
		type(_type), gcm_params(_gcm_params) {}

	virtual ~GCMCoupler() {}

	/** Called from within read_from_netcdf().  GCM-specific implementation
	sets up the contracts for this GCM - IceModel pair.  (Usually by calling
	through to GCM-specific virtual methods on the IceModel. */
	virtual void setup_contracts(
		IceModel &mod,
		NcFile &nc,
		std::string const &sheet_vname) = 0;

	/** @param sheets (OPTIONAL): IceSheet data structures w/ grids, etc. */
	virtual void read_from_netcdf(
		NcFile &nc, std::string const &vname,
		std::vector<std::string> const &sheet_names,
	    giss::MapDict<std::string, IceSheet> &sheets);

	/** Returns a unique rank number for each node in the parallel computation.
	Useful for debugging-type output. */
	int rank();

protected:
	/** @param time_s Time since start of simulation, in seconds
	Fields contained in the SMBMsg are INPUTS from the GCM.  They are
	therefore arranged according to gmc_inputs.  GCM inputs are converted
	into ice model inputs within IceModel::run_timestep(), which
	is called at the end of this method.
	@see gmc_inputs*/
	void call_ice_model(
		IceModel *model,
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
		giss::DynArray<SMBMsg> &sbuf);

};

}
