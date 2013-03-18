#include <glint2/IceModel_DISMAL.hpp>
#include <cstdio>

namespace glint2 {

/** Query all the ice models to figure out what fields they need */
void IceModel_DISMAL::get_required_fields(std::set<IceField> &fields)
{
	fields.insert(IceField::MASS_FLUX);
	fields.insert(IceField::ENERGY_FLUX);
}

/** @param index Index of each grid value.
@param vals The values themselves -- could be SMB, Energy, something else...
TODO: More params need to be added.  Time, return values, etc. */
void IceModel_DISMAL::run_timestep(
	blitz::Array<int,1> const &indices,
	std::map<IceField, blitz::Array<double,1>> const &vals2)
{
	printf("DISMAL: Run Timestep\n");
}


}
