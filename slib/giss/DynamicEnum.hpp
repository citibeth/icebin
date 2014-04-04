#pragma once

#include <string>
#include <vector>
#include <map>

namespace giss {


class DynamicEnum {
public:
	virtual ~DynamicEnum() {}

	virtual long size() const = 0;
	virtual int operator[](std::string const &name) const = 0;
	virtual std::string const &operator[](int ix) const = 0;
};


}
