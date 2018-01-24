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

#include <icebin/GCMRegridder.hpp>

namespace icebin {

BOOST_ENUM_VALUES( IceExch, int,
    (ICE)   (0)
    (EXCH)  (1)
)

class IceRegridder_L0 : public IceRegridder
{       // For ice model with level-value grid cells
public:
    /** The grid we use at the interpolation grid (exchange or ice) */
    IceExch interp_grid;

    /** Number of grid cells in the ice grid */
    size_t nI() const
        { return agridI.ndata(); }

    /** Number of grid cells in the interpolation grid */
    size_t nG() const
        { return interp_grid == IceExch::ICE ? nI() : aexgrid.ndata(); }

    IceRegridder_L0() : interp_grid(IceExch::EXCH) {}

public:
    // Implementations of virtual functions
    void GvEp(MakeDenseEigenT::AccumT &ret,
        blitz::Array<double,1> const *elevI) const;
    void GvI(MakeDenseEigenT::AccumT &ret,
        blitz::Array<double,1> const *elevI) const;
    void GvAp(MakeDenseEigenT::AccumT &ret,
        blitz::Array<double,1> const *elevI) const;
    void ncio(ibmisc::NcIO &ncio, std::string const &vname);
};


}   // namespace icebin
