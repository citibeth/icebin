#pragma once

#include <udunits2.h>
#include <iostream>

namespace giss {

class UTSystem;

class UTUnit
{
	friend class UTSystem;
	friend class CVConverter;
	friend std::ostream &operator<<(std::ostream &out, UTUnit const &unit);

	ut_unit *_self;
	bool _free_me;
	std::string _str;	// String used to instantiate this unit

	UTUnit(ut_unit *self, bool free_me, std::string const &str) :
		_self(self), _free_me(free_me), _str(str) {}

public:

	UTUnit() : _self(0), _free_me(false) {}

	~UTUnit() {
		if (_free_me && _self) ut_free(_self);
	}

	UTSystem get_system();

	std::string format(unsigned opts = UT_ASCII | UT_DEFINITION) const;

	std::string const &str() const
		{ return _str; }

	char const *c_str() const
		{ return _str.c_str(); }

	// ---------- Implement Move Semantics
#if 1
	UTUnit(UTUnit const &) = delete;
	UTUnit& operator=(UTUnit const &src) = delete;
#else
	UTUnit(UTUnit const &);
	UTUnit& operator=(UTUnit const &src);
#endif

	UTUnit(UTUnit &&src) {
		_self = src._self;
		_free_me = src._free_me;
		_str = std::move(src._str);
		src._self = 0;
	}

	UTUnit &operator=(UTUnit &&src) {
		_self = src._self;
		_free_me = src._free_me;
		_str = std::move(src._str);
		src._self = 0;
		return *this;
	}

};

inline std::ostream &operator<<(std::ostream &out, UTUnit const &unit)
	{ return out << unit.str(); }

#if 0
class UTException : public std::exception
{
	std::string _what;
public:
	UTException(std::string &&what) : _what(std::move(what)) {}
	UTException(std::string &what) : _what(what) {}
	UTException(boost::format &fmt) : _what(fmt.str()) {}
	const char *what() { return _what.c_str(); }
};
#endif


class UTSystem
{
	friend class UTUnit;

	ut_system *_self;
	bool _free_me;

	UTSystem(ut_system *self) : _self(self), _free_me(false) {}

public:
	UTSystem();

	UTSystem(std::string const &path);

	~UTSystem();

	UTUnit get_unit_by_name(std::string const &name) const;

	UTUnit get_unit_by_symbol(std::string const &symbol) const;

	UTUnit get_dimensionless_unit_one() const;

	UTUnit parse(std::string const &str, ut_encoding encoding = UT_ASCII) const;


	// ---------- Implement Move Semantics
	UTSystem(UTUnit const &) = delete;
	UTSystem& operator=(UTSystem const&) = delete;

	UTSystem(UTSystem &&src) {
		_self = src._self;
		src._self = 0;
	}

	UTSystem &operator=(UTSystem &&src) {
		_self = src._self;
		src._self = 0;
		return *this;
	}

};

class CVConverter
{
	cv_converter *_self;

public:
	CVConverter(UTUnit const &from, UTUnit const &to);

	~CVConverter()
		{ if (_self) cv_free(_self); }

	double convert(double const val) const
		{ return cv_convert_double(_self, val); }

	double *convert(double const *in, size_t count, double *out) const
		{ return cv_convert_doubles(_self, in, count, out); }


	// ---------- Implement Move Semantics
	CVConverter(UTUnit const &) = delete;
	CVConverter& operator=(CVConverter const&) = delete;

	CVConverter(CVConverter &&src) {
		_self = src._self;
		src._self = 0;
	}

	CVConverter &operator=(CVConverter &&src) {
		_self = src._self;
		src._self = 0;
		return *this;
	}

};

// =================================================
#if 0
inline UTSystem UTUnit::get_system()
	{ return UTSystem(ut_get_system(_self)); }

inline UTUnit& UTUnit::operator=(UTUnit const &src) {
	UTUnit dup(get_system().parse(src._str));
	_self = dup._self;
	_free_me = dup._free_me;
	_str = std::move(dup._str);
	dup._self = 0;
	return *this;
}

inline UTUnit::UTUnit(UTUnit const &src) {
	UTUnit dup(get_system().parse(src._str));
	_self = dup._self;
	_free_me = dup._free_me;
	_str = std::move(dup._str);
	dup._self = 0;
}
#endif


}
