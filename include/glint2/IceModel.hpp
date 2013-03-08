#pragma once

namespace glint2 {

/** The different things we can pass to an ice model. */
enum class IceField {
	MASS_FLUX,		// kg/(s m^2)
	ENERGY_FLUX		// W/m^2
};

class IceModel {

	BOOST_ENUM_VALUES( Type, int,
		(DISMAL)		(0)		// Demo Ice Sheet Model and LandIce
		(PISM)			(1)
		(ISSM)			(2)
	)

	/** Initialize any grid information, etc. from the IceSheet struct. */
	virtual void init(IceSheet *sheet);

	/** Query all the ice models to figure out what fields they need */
	virtual void get_required_fields(std::set<IceField> &fields) = 0;

	/** @param index Index of each grid value.
	@param vals The values themselves -- could be SMB, Energy, something else...
	TODO: More params need to be added.  Time, return values, etc. */
	virtual void run_timestep(
		blitz::Array<double,1> const &indices,
		std::map<IceField, blitz::Array<double,1>> const &vals2) = 0;

};

std::unique_ptr<IceModel> read_icemodel(NcFile &nc, std::string const &vname);

}
