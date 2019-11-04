#include <icebin/VarSet.hpp>
#include <ibmisc/netcdf.hpp>

namespace icebin {

int VarSet::add(
    std::string const &name,
    double default_value, std::string const &units,
    double mm, double bb,
    unsigned flags,
    std::string const &description)
{
    size_t ix = index.insert(name);

    VarMeta datum;
    datum.name = name;
    datum.default_value = default_value;
    datum.units = units;
    datum.flags = flags;
    datum.mm = mm;
    datum.bb = bb;
    datum.description = description;
    data.push_back(std::move(datum));
    return ix;
}

int VarSet::add(
    std::string const &name,
    double default_value, std::string const &units,
    unsigned flags,
    std::string const &description)
{
    size_t ix = index.insert(name);

    VarMeta datum;
    datum.name = name;
    datum.default_value = default_value;
    datum.units = units;
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
        ncvar.putAtt("units", var.units);
        ncvar.putAtt("description", var.description);
    }
}

}
