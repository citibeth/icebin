#pragma once

#include <string>
#include <vector>
#include <map>

namespace giss {

/** A two-way mapping between strings and ints (densely spaced), that
can be changed at runtime.  This serves as a virtual base class,
allowing algorithms to be written in terms of DynamicEnum (eg CouplingContract).

The string "unit" is special.  By convention, it is the LAST item in the enum.
Size can be queried with or without the "unit" string included. */
class DynamicEnum {
public:
	virtual ~DynamicEnum() {}

	/** Look up the int corresponding to a string. */
	virtual int operator[](std::string const &name) const = 0;
	/** Look up the string corresponding to an int. */
	virtual std::string const &operator[](int ix) const = 0;

	/** Number of strings in the enum, NOT COUNTING "unit" (if it exists) */
	virtual long size_withunit() const = 0;

	/** Number of strings in the enum, INCLUDING "unit" (if it exists) */
	virtual long size_nounit() const = 0;

	/** @return Index of the "unit" string.  Throws an exception if
	"unit" does not exist.
	NOTE: dynamic_enum[dynamc_enum.unit_ix()] == "unit" */
	virtual int unit_ix() const = 0;

};


}
