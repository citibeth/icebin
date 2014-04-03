#pragma once

#include <string>
#include <vector>
#include <map>

namespace giss {

/** A string-to-int mapping that is dynamic at runtime. */
class DynamicEnum {
	std::map<std::string, int> _name_to_ix;
	std::vector<std::string> _ix_to_name;

public:
	long size() { return _ix_to_name.size(); }

	/** Creates the enum. */
	void set_names(std::vector &&names);

	/** Given a string, returns its corresponding value.
	If the string is not in the enum, returns -1. */
	int operator[](std::string const &name)
	{
		auto ii = _name_to_ix.find(name);
		if (ii == _name_to_ix.end()) return -1;
		return *ii;
	}

	std::string const &operator[](int ix)
		{ return _ix_to_name[ix]; }


};

}
