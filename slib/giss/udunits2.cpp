#include <cstdio>
#include <cstring>
#include <exception>
#include <giss/udunits2.hpp>

namespace giss {

	std::string UTUnit::format(unsigned opts) const {
		char buf[1000];
		size_t len = ut_format(_self, buf, sizeof(buf), opts);
		return std::string(buf, len);
	}


	UTSystem::UTSystem() : _self(ut_new_system()), _free_me(true) {}

	UTSystem::UTSystem(std::string const &path) :
		_self(ut_read_xml(path == "" ? NULL : path.c_str())), _free_me(true) {}

	UTSystem::~UTSystem()
		{ if (_self && _free_me) ut_free_system(_self); }

	UTUnit UTSystem::get_unit_by_name(std::string const &name) const
	{
		UTUnit ret(ut_get_unit_by_name(_self, name.c_str()), false, name);
		if (ret._self) return ret;

		switch(ut_get_status()) {
			case UT_SUCCESS :
				fprintf(stderr,"UTUnit::get_unit_by_name(): Name %s doesn't map to a unit in the unit system.\n", name.c_str());
			break;
			case UT_BAD_ARG :
				fprintf(stderr, "UTUnit::get_unit_by_name(): UT_BAD_ARG, system or symbol is null.  This should not happen in C++ wrapper.\n");
			break;
			default :
				fprintf(stderr, "UTUnit::get_unit_by_name(): Unknown error\n");
			break;
		}
		throw std::exception();
	}

	UTUnit UTSystem::get_unit_by_symbol(std::string const &symbol) const
	{
		UTUnit ret(ut_get_unit_by_symbol(_self, symbol.c_str()), false, symbol);
		if (ret._self) return ret;

		switch(ut_get_status()) {
			case UT_SUCCESS :
				fprintf(stderr, "UTSystem::get_unit_by_system(): Symbol '%s' doesn't map to a unit in the unit system.\n", symbol.c_str());
			break;
			case UT_BAD_ARG :
				fprintf(stderr, "UTSystem::get_unit_by_system(): UT_BAD_ARG, system or symbol is null.  This should not happen in C++ wrapper.\n");
			break;
			default :
				fprintf(stderr, "UTSystem::get_unit_by_system(): Unknown error\n");
			break;
		}
		throw std::exception();
	}

	UTUnit UTSystem::get_dimensionless_unit_one() const
		{ return UTUnit(ut_get_dimensionless_unit_one(_self), true, "1"); }

	UTUnit UTSystem::parse(std::string const &str, ut_encoding encoding) const
	{
		char cstr[str.size()+1];
		strcpy(cstr, str.c_str());
		ut_trim(cstr, encoding);

		UTUnit ret(ut_parse(_self, cstr, encoding), true, str);
		if (ret._self) return ret;

		switch(ut_get_status()) {
			case UT_BAD_ARG :
				fprintf(stderr, "UTSystem::parse(): UT_BAD_ARG, system or str is null.  This should not happen in C++ wrapper.\n");
			break;
			case UT_SYNTAX :
				fprintf(stderr, "UTSystem::parse(): UT_SYNTAX error in '%s'\n", str.c_str());
			break;
			case UT_UNKNOWN :
				fprintf(stderr, "UTSystem::parse(): String '%s' contains an unknown identifier", str.c_str());
			break;
			case UT_OS :
				fprintf(stderr, "UTSystem::parse(): UT_OS\n");
			break;
			default :
				fprintf(stderr, "UTSystem::parse(): Unknown error\n");
			break;
		}
		throw std::exception();

	}



	CVConverter::CVConverter(UTUnit const &from, UTUnit const &to)
		: _self(ut_get_converter(from._self, to._self))
	{
		if (_self) return;

		switch(ut_get_status()) {
			case UT_BAD_ARG :
				fprintf(stderr, "CVConverter(%s -> %s): UT_BAD_ARG\n", from.c_str(), to.c_str()); break;
			case UT_NOT_SAME_SYSTEM :
				fprintf(stderr, "CVConverter(%s -> %s): UT_NOT_SAME_SYSTEM\n", from.c_str(), to.c_str()); break;
			case UT_MEANINGLESS :
				fprintf(stderr, "CVConverter(%s -> %s): UT_MEANINGLESS\n", from.c_str(), to.c_str()); break;
			case UT_OS :
				fprintf(stderr, "CVConverter(%s -> %s): UT_OS\n", from.c_str(), to.c_str()); break;
			default :
				fprintf(stderr, "CVConverter(%s -> %s): Unknown problem\n", from.c_str(), to.c_str()); break;
		}
		throw std::exception();

	}


}	// namespace giss

