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

class IceRegridder_L0 : public IceRegridder
{       // For ice model with level-value grid cells
public:
    /** Number of grid cells in the ice grid */
    size_t nI() const
        { return agridI.dim.sparse_extent(); }

    /** Number of grid cells in the interpolation grid */
    size_t nX() const
        { return aexgrid.sparse_extent(); }

    /** Number of grid cells in the Ice (I) or Exchange (X) grids.
    @param gridG Either 'I' or 'X' */
    size_t nG(char gridG) const;

public:
    // Implementations of virtual functions
    void GvEp(MakeDenseEigenT::AccumT &&ret,
        char gridG,    // Identity of G: 'I' (ice) or 'X' (exchange)
        blitz::Array<double,1> const *elevmaskI) const;
    void GvI(MakeDenseEigenT::AccumT &&ret,
        char gridG,    // Identity of G: 'I' (ice) or 'X' (exchange)
        blitz::Array<double,1> const *elevmaskI) const;
    void GvAp(MakeDenseEigenT::AccumT &&ret,
        char gridG,    // Identity of G: 'I' (ice) or 'X' (exchange)
        blitz::Array<double,1> const *elevmaskI) const;
    void ncio(ibmisc::NcIO &ncio, std::string const &vname);
};


}   // namespace icebin
