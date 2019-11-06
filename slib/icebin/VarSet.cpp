#include <icebin/VarSet.hpp>
#include <ibmisc/netcdf.hpp>

using namespace ibmisc;

namespace icebin {

double VarMeta::nc_factor(UTSystem const &ut_system) const
{
    // If no ncunits given, just use diagnostic units.
    if (this->ncunits.size() == 0) return 1.0;

    UTUnit unit0(ut_system.parse(units));
    UTUnit unit1(ut_system.parse(ncunits));

    // Obtain unit conversion factor
    try {
        CVConverter cv(unit0, unit1);
        return cv.convert(1.0);
    } catch(const std::exception &ex) {
        (*ibmisc_error)(-1,
            "Exception converting units for %s (%s -> %s)\n", name.c_str(), units.c_str(), ncunits.c_str());
    }
}


int VarSet::add(
    std::string const &name,
    double default_value,
    std::string const &units, std::string const &ncunits,
    double mm, double bb,
    unsigned flags,
    std::string const &description)
{
    size_t ix = index.insert(name);

    VarMeta datum;
    datum.name = name;
    datum.default_value = default_value;
    datum.units = units;
    datum.ncunits = ncunits;
    datum.flags = flags;
    datum.mm = mm;
    datum.bb = bb;
    datum.description = description;
    data.push_back(std::move(datum));
    return ix;
}

int VarSet::add(
    std::string const &name,
    double default_value,
    std::string const &units, std::string const &ncunits,
    unsigned flags,
    std::string const &description)
{
    size_t ix = index.insert(name);

    VarMeta datum;
    datum.name = name;
    datum.default_value = default_value;
    datum.units = units;
    datum.ncunits = ncunits;
    datum.flags = flags;
    datum.mm = 1.0;
    datum.bb = 0.0;
    datum.description = description;
    data.push_back(std::move(datum));
    return ix;
}

netCDF::NcVar VarSet::ncdefine(ibmisc::NcIO &ncio,
    std::vector<netCDF::NcDim> const &dims,
    std::string vname_base) const
{
    for (size_t i=0; i<size(); ++i) {
        VarMeta const &var((*this)[i]);
        auto ncvar(get_or_add_var(ncio, vname_base+var.name, "double", dims));
        ncvar.putAtt("units", var.final_nc_units());
        ncvar.putAtt("description", var.description);
    }
}



}
