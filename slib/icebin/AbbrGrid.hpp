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

#include <ibmisc/enum.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/iter.hpp>
#include <ibmisc/Proj2.hpp>
#include <ibmisc/indexing.hpp>

#include <icebin/error.hpp>

#ifdef BUILD_MODELE
#include <icebin/modele/hntr.hpp>
#endif

namespace icebin {

class AbbrGrid {
    Grid::Type type;
    std::shared_ptr<GridExtra> extra;
    Grid::Coordinates coordinates;
    Grid::Parameterization parameterization;
    std::string name;
    std::string sproj;


    // Items here will correspond to cells (for L0) or vertices (for L1)
    SparseSet<long,int> dim;
    std::array<blitz::Array<int,1>,3> ijk;
    blitz::Array<double,1> native_area;    // dense indexing
    blitz::Array<double,1> proj_area;
    blitz::Array<double,1> centroid;

    size_t ndata() { return dim.sparse_extent(); }

    virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname, bool rw_full=true);

    AbbrGrid(Grid &g);

};

AbbrGrid::AbbrGrid(Grid &g)
{
    type = g.type;
    extra = g.extra;
    coordinates = g.coordinates;
    parameterization = g.parameterization;
    name = g.name;
    sproj = g.sproj;

    // Allocate
    nd = g.ndata();
    for (int i=0; i<3; ++i) ijk[i].reference(blitz::Array<int,1>(nd));
    native_area.reference(blitz::Array<double,1>(nd));
    proj_area.reference(blitz::Array<double,1>(nd));
    centroid.reference(blitz::Array<Point,1>(nd));

    // Copy info into AbbrGrid
    int id=0;
    for (auto cell=g.cells.begin(); cell != g.cells.end(); ++cell) {
        dim.add_dense(cell->index);
        ijk[0](id) = cell->i;
        ijk[1](id) = cell->j;
        ijk[2](id) = cell->k;
        native_area(id) = cell->native_area;
        proj_area(id) = cell->proj_area;
        centroid(id) = cell->centroid;
    }
}

}    // namespace

