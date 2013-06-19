#pragma once

#include <glint2/IceModel.hpp>

namespace glint2 {

class IceModel_Decode : public IceModel {
public:
	// Dimensionality Ice Model's vector space
	int const ndata;

	IceModel_Decode(Grid const &grid) : ndata(grid.ndata()) {}
	IceModel_Decode(int _ndata) : ndata(_ndata) {}

	/** @param index Index of each grid value.
	@param vals The values themselves -- could be SMB, Energy, something else...
	TODO: More params need to be added.  Time, return values, etc.
	@param itime Some kind of representation of the current GCM timestep.
	Helps with debugging. */
	virtual void run_timestep(long itime,
		blitz::Array<int,1> const &indices,
		std::map<IceField, blitz::Array<double,1>> const &vals2);

	/** Runs a timestep after fields have been decoded.  This is what
	one will normally want to override, unless you wish to decode
	yourself. */
	virtual void run_decoded(long itime,
		std::map<IceField, blitz::Array<double,1>> const &vals2) = 0;

};

}
