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

class ExchangeGrid {
    // Sparse indexing needed by IceRegridder::init()
    std::vector<int> indices;    // Length*2: (ixB, ixA)
    std::vector<double> overlaps;

public:
    ExchangeGrid() {}

    /** Only works for Grid objects resulting from the overlap program. */
    explicit ExchangeGrid(Grid const &g);

    void reserve(size_t n)
    {
        indices.reserve(n*2);
        overlaps.reserve(n);
    }

    void add(std::array<int,2> const &index, double _area)
    {
        indices.push_back(index[0]);
        indices.push_back(index[1]);
        overlaps.push_back(_area);
    }

    int dense_extent() const 
        { return overlaps.size(); }

    long sparse_extent() const
        { return overlaps.size(); }

    /** Exchange gridcells are numbered in order from 0.
    Therefore, dense and sparse indexing are equivalent. */
    long to_sparse(int id) const
        { return id; }

    int ijk(int id, int index) const
        { return indices[id*2 + index]; }
    double native_area(int id) const
        { return overlaps[id]; }

    void ncio(ibmisc::NcIO &ncio, std::string const &vname);

    /** NOTE: This will result in ExchangeGrid cells being renumbered,
    resulting in different numbering schemes for different processors.
    That is not a problem because matrices based on this grid are only
    used temporarily; and this dimension is ultimately multiplied away
    before being shared between processors or in time. */
    void filter_cellsB(std::function<bool(long)> const &keep_B_fn);

};



struct AbbrGrid {
    ibmisc::clonable_unique_ptr<GridSpec> spec;
    GridCoordinates coordinates;
    GridParameterization parameterization;
    ibmisc::Indexing indexing;
    std::string name;
    std::string sproj;    // Set for ice grids, not GCM grid


    // Items here will correspond to cells (for L0) or vertices (for L1)
    spsparse::SparseSet<long,int> dim;    // cell->index = sparse
    blitz::Array<int,2> ijk;    // ijk(index, ijk)  (OPTIONAL; and second dimension varies in size)
    blitz::Array<double,1> native_area;    // dense indexing
    // blitz::Array<double,1> proj_area;

    // Only set if coordinates == GridCoordinates::XY
    blitz::Array<double,2> centroid_xy;    // centroid(index, xy)

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


    // ===============================================================
    // We only need to define these because blitz::Array does not follow STL conventions
    // and is not movable.
    void operator=(AbbrGrid &&other);

    AbbrGrid(AbbrGrid &&other);

    void operator=(AbbrGrid const &other);

    AbbrGrid(AbbrGrid const &other);

    // ===============================================================
};




}    // namespace

