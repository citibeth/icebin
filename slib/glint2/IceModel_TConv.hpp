#pragma once

#include <glint2/IceModel_Decode.hpp>

namespace glint2 {

/** Produces a Surface T field for ice models that need it,
from GCM-supplied Mass and Energy fields */
class IceModel_TConv : public IceModel_Decode {
	std::unique_ptr<IceModel_Decode> model;

	double const LHM;	// latent heat of melt at 0 C (334590 J/kg)
	double const SHI;	// heat capacity of pure ice (at 0 C) (2060 J/kg C)

public:
	IceModel_TConv(std::unique_ptr<IceModel_Decode> &&_model,
		double LHM, double SHI);

 	/** Query all the ice models to figure out what fields they need */
	void get_required_fields(std::set<IceField> &fields);

	/** @param index Index of each grid value.
	@param vals The values themselves -- could be SMB, Energy, something else...
	TODO: More params need to be added.  Time, return values, etc.
	@param itime Some kind of representation of the current GCM timestep.
	Helps with debugging. */
	void run_decoded(long itime,
		std::map<IceField, blitz::Array<double,1>> const &vals2);

	void read_from_netcdf(NcFile &nc, std::string const &vname);



};


}