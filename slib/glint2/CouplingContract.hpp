#pragma once

#include <giss/cfnames.hpp>
#include <giss/DynamicEnum.hpp>
#include <giss/VarMetaData.hpp>
#include <ostream>

namespace glint2 {

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
	units(_units),
	description(_description)
	{}


	CoupledField(giss::CFName const *cf) :
		name(cf->id), description(cf->description), units(cf->canonical_units)
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
	CouplingContract() : _size_nounit(0), _unit_ix(-1) {}

	auto begin() const
		{ return _ix_to_field.begin(); }
	auto end() const
		{ return _ix_to_field.end(); }


	void add_field(CoupledField const &cf)
	{
		int index = _ix_to_field.size();
		_name_to_ix.insert(std::make_pair(cf.name, index));
		_ix_to_field.push_back(cf);
		if (cf.name == "unit") {
			_unit_ix = index;
		} else {
			_size_nounit += 1;
		}
	}

	void add_field(giss::CFName const *cf)
		{ add_field(CoupledField(cf->id, cf->canonical_units, cf->description)); }

	void add_field(std::string const &name, std::string const &units,
		std::string const &description)
	{ add_field(CoupledField(name, units, description)); }


	void add_cfname(std::string const &name, std::string const &units = "")
		{ add_field(giss::get_cfname(name, units)); }


	long size_withunit() const { return _ix_to_field.size(); }
	long size_nounit() const { return _size_nounit; }

	int unit_ix() const { return _unit_ix; }

	int operator[](std::string const &name) const {
		auto ii = _name_to_ix.find(name);
		if (ii == _name_to_ix.end()) return -1;
		return ii->second;
	}

	std::string const &operator[](int ix) const
		{ return _ix_to_field[ix].name; }

	std::ostream &operator<<(std::ostream &out) const;
 
	friend std::ostream &operator<<(std::ostream &out, CouplingContract const &con);

};

extern std::ostream &operator<<(std::ostream &out, CouplingContract const &con);


}
