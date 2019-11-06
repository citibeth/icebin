/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <string>
#include <ibmisc/IndexSet.hpp>
#include <ibmisc/VarTransformer.hpp>
#include <ibmisc/udunits2.hpp>

namespace netCDF {
    class NcVar;
    class NcDim;
}
namespace ibmisc {
    class NcIO;
}

namespace icebin {


struct VarMeta {
    /** The "short" name of the variable */
    std::string name;
    /** The units of the variable, in UDUNITS format.
    NOTE: This is the units BEFORE rescaling and sending to the GCM. */
    std::string units;          //!< UDUnits-compatible string
    /** The units to convert to when writing NetCDF (BEFORE rescaling) */
    std::string ncunits;
    /** The flags the variable resides on; see contracts.hpp for values */
    unsigned flags;         //!< Allows arbitrary subsets
    /** A textual description of the variable, also called the "long name" */
    std::string description;
    /** Rescaling factors... used in apply_gcm_inputs().  result =
        mm*source + bb */
    double mm = 1;
    double bb = 0;

    double default_value;

    /** Units to use in NetCDF file */
    std::string const &final_nc_units() const
        { return ncunits.size() == 0 ? units : ncunits; }

    /** Computes the conversion factor from units -> ncunits */
    double nc_factor(ibmisc::UTSystem const &ut_system) const;

};

class VarSet
{
public:

    ibmisc::IndexSet<std::string> index;    // Densely ordered set of constant names
    std::vector<VarMeta> data;  // Meta-data and value data

    std::vector<std::string> const &keys() const
        { return index.keys(); }

    int add(
        std::string const &name,
        double default_value,
        std::string const &units, std::string const &ncunits,
        double mm, double bb,
        unsigned flags = 0,
        std::string const &description = "<no description>");

    int add(
        std::string const &name,
        double default_value,
        std::string const &units, std::string const &ncunits,
        unsigned flags = 0,
        std::string const &description = "<no description>");

    size_t size() const { return index.size(); }

    VarMeta const &operator[](size_t ix) const
        { return data[ix]; }

    VarMeta const &at(std::string const &name) const
        { return data[index.at(name)]; }

    /** Defines each variable in the DimSet, according to a set of dimensions. */
    netCDF::NcVar ncdefine(ibmisc::NcIO &ncio,
        std::vector<netCDF::NcDim> const &dims,
        std::string vname_base) const;

};


}   // Namespace

inline std::ostream &operator<<(std::ostream &out, icebin::VarMeta const &cf)
    { return out << "(" << cf.name << ": [" << cf.units << "] flags:" << cf.flags << ")"; } 

inline std::ostream &operator<<(std::ostream &out, icebin::VarSet const &vars)
{
    for (auto ii=vars.data.begin(); ii != vars.data.end(); ++ii) {
        std::cout << *ii << std::endl;
    }
    return out;
}

