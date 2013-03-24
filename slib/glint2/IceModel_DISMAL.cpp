#include <glint2/IceModel_DISMAL.hpp>
#include <cstdio>
#include <cmath>

namespace glint2 {

/** Query all the ice models to figure out what fields they need */
void IceModel_DISMAL::get_required_fields(std::set<IceField> &fields)
{
	fields.insert(IceField::MASS_FLUX);
	fields.insert(IceField::ENERGY_FLUX);
}


static double const nan = std::numeric_limits<double>::quiet_NaN();

blitz::Array<double,2> IceModel_DISMAL::decode(
	blitz::Array<int,1> const &indices,
	blitz::Array<double,1> const &vals)
{
	blitz::Array<double,2> ret(ny, nx);
	ret = nan;

	// Make a 1-D alias of this array, and copy into it
	int extent = ny*nx;
	blitz::Array<double,1> ret_1d(ret.data(),
		blitz::shape(extent), blitz::neverDeleteData);

	int n = indices.size();
	for (int i=0; i < n; ++i) {
		int ix = indices(i);
		// Do our own bounds checking!
		if (ix < 0 || ix >= extent) {
			fprintf(stderr, "IceModel_DISMAL: index %d out of range [0, %d)\n", ix, extent);
			throw std::exception();
		}

		// Sanity check for NaN coming through
		if (std::isnan(vals(i))) {
			fprintf(stderr, "IceModel_DISMAL::decode: vals[%d] (index=%d) is NaN!\n", i, ix);
			throw std::exception();
		}
		double &oval = ret_1d(ix);
		if (std::isnan(oval)) oval = vals(i);
		else oval += vals(i);
	}

	return ret;
}


/** @param index Index of each grid value.
@param vals The values themselves -- could be SMB, Energy, something else...
TODO: More params need to be added.  Time, return values, etc. */
void IceModel_DISMAL::run_timestep(
	blitz::Array<int,1> const &indices,
	std::map<IceField, blitz::Array<double,1>> const &vals2)
{
	printf("DISMAL: Run Timestep\n");
	auto mass(decode(indices, vals2.find(IceField::MASS_FLUX)->second));
	auto energy(decode(indices, vals2.find(IceField::ENERGY_FLUX)->second));

	NcFile ncout("dismal.nc", NcFile::Replace);

	std::vector<boost::function<void ()>> fns;
	NcDim *nx_dim = ncout.add_dim("nx", nx);
	NcDim *ny_dim = ncout.add_dim("ny", ny);

	fns.push_back(giss::netcdf_define(ncout,
		"mass", mass, {ny_dim, nx_dim}));

	fns.push_back(giss::netcdf_define(ncout,
		"energy", energy, {ny_dim, nx_dim}));

	for (auto ii = fns.begin(); ii != fns.end(); ++ii) (*ii)();

	ncout.close();


}


}
