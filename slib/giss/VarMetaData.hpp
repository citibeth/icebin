#pragma once

#include <string>

namespace giss {

/** Virtual base class providing meta-data about variables in a system.
The idea is to decouple this meta-data from the actual data storage format,
which might vary across frameworks (eg. ModelE, PetSC, blitz++, etc).
This meta-data can be used when writing out stuff to NetCDF. */
class VarMetaData {
public:
	virtual ~VarMetaData() {}

	/** The "short" name of the variable */
	virtual std::string const &get_name() const = 0;

	/** The units of the variable, in UDUNITS format. */
	virtual std::string const &get_units() const { return ""; }

	/** A textual description of the variable, also called the "long name" */
	virtual std::string const &get_description() const { return ""; }
};

}
