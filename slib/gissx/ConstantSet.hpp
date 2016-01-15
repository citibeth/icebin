#pragma once

#include <giss/CouplingContract.hpp>
#include <giss/udunits2.hpp>

class NcFile;

namespace giss {

class ConstantSet {
	CouplingContract fields;
	std::vector<double> vals;

public:
	UTSystem const * const ut_system;	//!< Unit system to use for conversions

	ConstantSet(UTSystem const *_ut_system) : ut_system(_ut_system) {}

	int add_field(
		std::string const &name, std::string const &units,
		std::string const &description);

	/** @return Index of the new constant. */
	int set(
		std::string const &name,
		double val,
		std::string const &units,
		std::string const &description);

	std::string const &units(std::string const &name)
		{ return fields.field(name).units; }

	/** Copy a constant's value and units to another constant
	@return Index of the new constant. */
	int copy(std::string const &dst_name,
		std::string const &src_name,
		std::string const &description);

#if 0
	/** Transfers a constant from a different constant set into
	this one, converting between units if needed. */
	double set(std::string const &name,
		ConstantSet const &src,
		std::string const &src_name);
#endif

	double get_as(std::string const &name,
		UTUnit const &units) const;

	double get_as(std::string const &name,
		std::string const &sunits) const;


	double const &operator[](std::string const &name) const
		{ return vals[fields.index(name)]; }

	double const &operator[](int ix) const
		{ return vals[ix]; }

	double &operator[](std::string const &name)
		{ return vals[fields.index(name)]; }

	double &operator[](int ix)
		{ return vals[ix]; }

	std::vector<CoupledField>::const_iterator begin()
		{ return fields.begin(); }
	std::vector<CoupledField>::const_iterator end()
		{ return fields.end(); }

	void netcdf_define(NcFile &nc, std::string const &vname);
	void read_from_netcdf(NcFile &nc, std::string const &vname);
};

}
