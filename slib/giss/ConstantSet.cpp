#include <cstring>
#include <giss/ConstantSet.hpp>
#include <giss/udunits2.hpp>
#include <giss/ncutil.hpp>
#include <giss/string.hpp>
#include <giss/exit.hpp>

namespace giss {


int ConstantSet::add_field(
	std::string const &name, std::string const &units,
	std::string const &description)
{
	int ix = fields.add_field(name, units, 0, description);
	if (ix >= vals.size()) vals.resize(ix+1);
	return ix;
}

int ConstantSet::set(
	std::string const &name,
	double val,
	std::string const &units,
	std::string const &description)
{
	int ix = add_field(name, units, description);
	vals[ix] = val;
	return ix;
}

/** Copy a constant's value and units to another constant
@return Index of the new constant. */
int ConstantSet::copy(std::string const &dst_name,
	std::string const &src_name,
	std::string const &description)
{
	int src_ix = fields.index(src_name);
	CoupledField const &src_field(fields.field(src_ix));
	
	int dst_ix = add_field(dst_name, src_field.units, description);
	vals[dst_ix] = vals[src_ix];

	return dst_ix;
}


double ConstantSet::get_as(std::string const &name,
	UTUnit const &units) const
{
	int src_ix = fields.index(name);
	std::string const &susrc(fields.field(src_ix).units);
	UTUnit usrc(ut_system->parse(susrc));

	try {
		CVConverter cv(usrc, units);
		double ret = cv.convert((*this)[src_ix]);
printf("ConstantSet: Converting %s: %f %s --> %f %s\n", name.c_str(), (*this)[src_ix], usrc.c_str(), ret, units.c_str());
		return ret;
	} catch(const std::exception &ex) {
		fprintf(stderr, "Exception in ConstantSet::get_as(%s, %s)\n", name.c_str(), units.c_str());
		giss::exit(1);
	}
	
}

double ConstantSet::get_as(std::string const &name,
	std::string const &sunits) const
{
	UTUnit units(ut_system->parse(sunits));
	return get_as(name, units);
}


#if 0
/** Transfers a constant from a different constant set into
this one, converting between units if needed. */
double ConstantSet::set(std::string const &name,
	ConstantSet const &src,
	std::string const &src_name)
{
	int src_ix = src.fields[src_name];
	std::string const &susrc(src.fields.field(src_ix).units);
	UTUnit usrc(ut_system->parse(susrc));

	int dst_ix = fields[name];
	std::string const &sudst(fields.field(dst_ix).units);
	UTUnit udst(ut_system->parse(sudst));

	try {
		CVConverter cv(usrc, udst);
		double ret = cv.convert(src[src_ix]);
		(*this)[dst_ix] = ret;
printf("ConstantSet: Converting %f %d --> %f %d\n", src[src_ix], usrc.c_str(), ret, udst.c_str());
		return ret;
	} catch(const std::exception &ex) {
		fprintf(stderr, "Exception in ConstantSet::set(%s, ConstantSet, %s)\n", name.c_str(), src_name.c_str());
		giss::exit(1);
	}
}
#endif


void ConstantSet::netcdf_define(NcFile &nc, std::string const &vname)
{
	// Create the variable to store our constants
	NcVar *ncvar = giss::get_var_safe(nc, vname.c_str(), false);
	if (!ncvar) {
		auto oneDim = giss::get_or_add_dim(nc, "one", 1);
		ncvar = nc.add_var(vname.c_str(), ncDouble, oneDim);
	}

	// Store the constants as attributes
	for (int i=0; i<fields.size_withunit(); ++i) {
		CoupledField const &field(fields.field(i));
		ncvar->add_att(field.name.c_str(), vals[i]);
		ncvar->add_att((field.name + "_units").c_str(), field.units.c_str());
		ncvar->add_att((field.name + "_description").c_str(), field.description.c_str());
	}
}

void ConstantSet::read_from_netcdf(NcFile &nc, std::string const &vname)
{
	NcVar *ncvar = giss::get_var_safe(nc, vname.c_str(), true);
	if (!ncvar) {
		fprintf(stderr, "ConstantSet::read_from_netcdf() cannot find variable %s\n", vname.c_str());
		giss::exit(1);
	}
	int n = ncvar->num_atts();

	// Read through the attributes, getting units and names separately
	std::map<std::string, std::string> units;
	std::map<std::string, double> consts;
	std::map<std::string, std::string> descriptions;
	for (int i=0; i<n; ++i) {
		auto att = giss::get_att(ncvar, i);
		std::string att_name(att->name());
		if (giss::ends_with(att_name, "_description")) {
			descriptions.insert(std::make_pair(
				att_name.substr(0, att_name.size() - std::strlen("_description")),
				std::string(att->as_string(0))));
		} else if (giss::ends_with(att_name, "_units")) {
			units.insert(std::make_pair(
				att_name.substr(0, att_name.size() - std::strlen("_units")),
				std::string(att->as_string(0))));
		} else {
			consts.insert(std::make_pair(
				att_name, att->as_double(0)));
		}
	}

	// Now go through them again, matching up constants and units
	for (auto ii = consts.begin(); ii != consts.end(); ++ii) {
		std::string const &name = ii->first;
		double const val = ii->second;

		auto ui = units.find(name);
		if (ui == units.end()) {
			fprintf(stderr, "Could not find _units attribute for %s\n", name.c_str());
		}

		auto di = descriptions.find(name);
		if (di == descriptions.end()) {
			fprintf(stderr, "Could not find _description attribute for %s\n", name.c_str());
		}

		std::string const &u = units.find(name)->second;
		std::string const &d = descriptions.find(name)->second;

		set(name, val, u, d);
	}
}

}
