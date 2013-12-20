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
#include <glint2/IceModel.hpp>
#include <boost/filesystem.hpp>

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
protected:

	/** The MPI Communicator used by the GCM */
	MPI_Comm comm;

	/** Rank of the root process, which will receive data on the ice grid */
	int root;

public:
	// Only needed by root MPI node in MPI version
	giss::MapDict<int,IceModel> models;

	GCMCoupler(MPI_Comm _comm, int _root) :
		comm(_comm), root(_root) {}


	/** Query all the ice models to figure out what fields they need */
	std::set<IceField> get_required_fields();

	/** @param sheets (OPTIONAL): IceSheet data structures w/ grids, etc. */
	virtual void read_from_netcdf(
		boost::filesystem::path const &glint2_config_dir,
		NcFile &nc, std::string const &vname,
		std::vector<std::string> const &sheet_names,
	    giss::MapDict<std::string, IceSheet> const &sheets);

	/** Returns a unique rank number for each node in the parallel computation.
	Useful for debugging-type output. */
	virtual int rank();

protected:
	/** @param time_s Time since start of simulation, in seconds */
	void call_ice_model(double time_s,
		giss::DynArray<SMBMsg> &rbuf,
		std::vector<IceField> const &fields,
		SMBMsg *begin, SMBMsg *end);


public:
	/** @param sbuf the (filled) array of ice grid values for this MPI node.
	@param time_s Time (seconds) since the start of the GCM run.
	*/
	void couple_to_ice(double time_s,
		std::vector<IceField> const &fields,
		giss::DynArray<SMBMsg> &sbuf);

};

}
