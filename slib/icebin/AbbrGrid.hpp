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

#include <vector>
#include <unordered_map>
#include <functional>

#include <ibmisc/memory.hpp>
#include <ibmisc/enum.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/iter.hpp>
#include <ibmisc/Proj2.hpp>
#include <ibmisc/indexing.hpp>
#include <spsparse/SparseSet.hpp>

#include <icebin/error.hpp>
#include <icebin/GridSpec.hpp>

namespace icebin {

class Grid;

struct AbbrGrid {
    ibmisc::clonable_unique_ptr<GridSpec> spec;
    GridCoordinates coordinates;
    GridParameterization parameterization;
    ibmisc::Indexing indexing;
    std::string name;
    std::string sproj;    // Set for ice grids, not GCM grid


    // Items here will correspond to cells (for L0) or vertices (for L1)
    spsparse::SparseSet<long,int> dim;    // cell->index = sparse
    blitz::Array<int,2> ijk;    // ijk(index, ijk)
    blitz::Array<double,1> native_area;    // dense indexing
    // blitz::Array<double,1> proj_area;

    // Only set if coordinates == GridCoordinates::XY
    blitz::Array<double,2> centroid_xy;    // centroid(index, xy)

    size_t ndata() const { return dim.sparse_extent(); }

    virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname);

    AbbrGrid() {}
    explicit AbbrGrid(Grid const &g);


    AbbrGrid(
        std::unique_ptr<GridSpec> &&_spec,
        GridCoordinates _coordinates,
        GridParameterization _parameterization,
        ibmisc::Indexing _indexing,
        std::string const &_name,
        std::string const &_sproj,
        spsparse::SparseSet<long,int> &&_dim,
        blitz::Array<int,2> const &_ijk,
        blitz::Array<double,1> const &_native_area,
        blitz::Array<double,2> const &_centroid_xy);


    void clear();

    void filter_cells(std::function<bool(long)> const &keep_fn);



#if 0
    void operator=(AbbrGrid const &other)
    {
        type = other.type;
        spec = other.spec->clone();
        coordinates = other.coordinates;
        parameterization = other.parameterization;
        indexing = other.indexing;
        name = other.name;
        sproj = other.sproj;

        dim = other.dim;
        ijk.reference(other.ijk);
        native_area.reference(other.native_area);
        centroid.reference(other.centroid);
    }
#endif

};




}    // namespace

