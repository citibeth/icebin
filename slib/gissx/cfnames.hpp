#pragma once

#include <string>
#include <map>

namespace giss {

	/** Common meta-data for variables.  Can be written to NetCDF. */
	struct CFName {

		/** CF-compliant standard name.
		@see http://cf-pcmdi.llnl.gov/documents/cf-standard-names */
		std::string const id;
		std::string const description;
		/** UDUnits-compatible string */
		std::string const canonical_units;
		std::string const grib;
		std::string amip;
	};

	namespace cf {
		extern std::map<std::string, giss::CFName> const &standard_name_table();
	}	// namespace cf

	/** Retrieves a CFName record by name.  Throws an exception if the name is not
	found, or if the units of the CFName record is not equal to expected_units. */
	extern CFName const *get_cfname(
		std::string const &_id,
		std::string const &expected_units = "");

}	// namespace giss
