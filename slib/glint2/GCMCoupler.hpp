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
#include <giss/cfnames.hpp>

namespace glint2 {


struct CoupledField {
	std::string name;
	std::string description;
	std::string units;			//!< UDUnits-compatible string

	CoupledField(std::string const &_name, std::string const &_description,
		std::string const _&units) : name(_name), description(_description), units(_units) {}


	CoupledField(CFName *cfname) :
		name(cf->id), description(cf->description), units(cf->canonical_units)
	{}

	std::ostream operator<<(std::ostream &out)
		{ return out << "(" << name << ": " << units << ")"; } 

};

struct CouplingContract : public DynamicEnum
{
	std::vector<CoupledField> _ix_to_field;
	std::map<std::string, int> _name_to_ix;
public:

	void push_back(CoupledField const &cf)
	{
		_name_to_ix.insert(std::make_pair(cf.id, _ix_to_field.size()));
		_ix_to_field.push_back(cf);
	}

	long size() { return _ix_to_field.size(); }

	int operator[](std::string const &name) {
		auto ii = _name_to_ix.find(name);
		if (ii == _name_to_ix.end()) return -1;
		return *ii;
	}

	std::string const &operator[](int ix)
		{ return _ix_to_field[ix].id; }

	std::ostream operator<<(std::ostream &out);
};


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

	/** Items held in GCMCoupler on a per-ice-sheet basis. */
	struct PerSheet {
		std::unique_ptr<IceModel> model;

		/** Ordered specification of the variables (w/ units)
		to be passed from GLINT2 to the ice model */
		CouplingContract con_gcm_to_ice;

		/** Variable transformations to prepare GCM
		output for consumption by the ice model.
		This transformation is in addition to regridding. */
		VarTransformer vt_gcm_to_ice;

		/** Ordered specification of the variables (w/ units)
		to be passed from the ice model to GLINT2 */
		CouplingContract con_ice_to_gcm;

		/** Variable transformations to prepare ice model 
		output for consumption by the GCM.
		This transformation is in addition to regridding. */
		VarTransformer vt_ice_to_gcm;
	};

	IceModel::GCMParams const gcm_params;

	// Only needed by root MPI node in MPI version
	std::map<int, PerIce> per_sheet;

	/** Fields we receive from the GCM */
	CouplingContract gcm_inputs;

	GCMCoupler(Type _type, IceModel::GCMParams const &_gcm_params) :
		type(_type), gcm_params(_gcm_params) {}

	/** @param sheets (OPTIONAL): IceSheet data structures w/ grids, etc. */
	virtual void read_from_netcdf(
		NcFile &nc, std::string const &vname,
		std::vector<std::string> const &sheet_names,
	    giss::MapDict<std::string, IceSheet> &sheets);

	/** Returns a unique rank number for each node in the parallel computation.
	Useful for debugging-type output. */
	int rank();

protected:
	/** @param time_s Time since start of simulation, in seconds */
	void call_ice_model(
		IceModel *model,
		double time_s,
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
