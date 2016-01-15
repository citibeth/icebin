#pragma once

#include <string>

namespace icebin {

/** Virtual base class providing meta-data about variables in a system.
The idea is to decouple this meta-data from the actual data storage format,
which might vary across frameworks (eg. ModelE, PetSC, blitz++, etc).
This meta-data can be used when writing out stuff to NetCDF.

This class is not really general (for all meta-data in CF
conventions); more of a throwaway class for this application.
*/
class VarMetaData {
	std::string _name;
	std::string _units;			//!< UDUnits-compatible string
	unsigned _flags;			//!< Allows arbitrary subsets
	std::string _description;
	double _default_value;

public:

	VarMetaData(std::string const &name,
		double default_value,
		std::string const &units,
		unsigned flags,
		std::string const &description)
	: _name(name), _default_value(default_value),
	_units(units),
	_flags(flags),
	_description(description)
	{}

	virtual ~VarMetaData() {}


	/** The "short" name of the variable */
	std::string const &name() const { return _name; }

	/** The units of the variable, in UDUNITS format. */
	std::string const &units() const { return _units; }

	/** The flags the variable resides on. */
	unsigned flags() const { return _flags; }

	/** A textual description of the variable, also called the "long name" */
	std::string const &description() const { return _description; }

	double default_value() const { return _default_value; }
};

}
