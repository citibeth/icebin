#pragma once

#include <giss/cfnames.hpp>
#include <giss/DynamicEnum.hpp>
#include <giss/VarMetaData.hpp>
#include <ostream>

namespace giss {

struct CoupledField : public giss::VarMetaData {
	std::string name;
	std::string units;			//!< UDUnits-compatible string
	std::string description;

	// Implementing giss::VarMetaData
	std::string const &get_name() const { return name; }
	std::string const &get_units() const { return units; }
	std::string const &get_description() const { return description; }

	CoupledField(std::string const &_name,
		std::string const &_units,
		std::string const &_description)
	: name(_name),
	units(std::move(_units)),
	description(_description)
	{}

};

inline std::ostream &operator<<(std::ostream &out, CoupledField const &cf)
	{ return out << "(" << cf.name << ": " << cf.units << ")"; } 




struct CouplingContract : public giss::DynamicEnum
{
	std::vector<CoupledField> _ix_to_field;
	std::map<std::string, int> _name_to_ix;
	long _size_nounit;		//!< Num names, not including "unit"
	int _unit_ix;

public:
	CouplingContract() : _size_nounit(0), _unit_ix(1) {}

	auto begin() const -> decltype(_ix_to_field.begin())
		{ return _ix_to_field.begin(); }
	auto end() const -> decltype(_ix_to_field.end())
		{ return _ix_to_field.end(); }

	int add_field(CoupledField &&cf);

	int add_field(CoupledField const &cf) {
		CoupledField dup(cf);
		return add_field(std::move(cf));
	}

	int add_field(giss::CFName const *cf) {
		return add_field(CoupledField(
			cf->id, cf->canonical_units, cf->description));
	}

	int add_field(std::string const &name, std::string const &units,
		std::string const &description)
	{ return add_field(CoupledField(name, units, description)); }


	int add_cfname(std::string const &name, std::string const &units = "")
		{ return add_field(giss::get_cfname(name, units)); }


	long size_withunit() const { return _ix_to_field.size(); }
	long size() const { return size_withunit(); }
	long size_nounit() const { return _size_nounit; }

	int unit_ix() const { return _unit_ix; }

	int operator[](std::string const &name) const {
		auto ii = _name_to_ix.find(name);
		if (ii == _name_to_ix.end()) return -1;
		return ii->second;
	}

	std::string const &operator[](int ix) const
		{ return _ix_to_field[ix].name; }

	CoupledField const &field(int ix) const
		{ return _ix_to_field[ix]; }

	CoupledField const &field(std::string const &name) const
		{ return field((*this)[name]); }

	std::ostream &operator<<(std::ostream &out) const;
 
	friend std::ostream &operator<<(std::ostream &out, CouplingContract const &con);

};

extern std::ostream &operator<<(std::ostream &out, CouplingContract const &con);


}
