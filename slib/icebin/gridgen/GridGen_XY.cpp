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

#include <functional>
#include <cmath>
#include <ibmisc/netcdf.hpp>

#include <icebin/gridgen/GridGen_XY.hpp>
#include <icebin/gridgen/gridutil.hpp>

using namespace netCDF;
using namespace ibmisc;

namespace icebin {

GridSpec_XY GridSpec_XY::make_with_boundaries(
    double x0, double x1, double dx,
    double y0, double y1, double dy)
{
    GridGen_XY spec;

    // Set up x coordinates
    int nx = (int)(.5 + (x1 - x0) / dx);    // Round to nearest integer
    double nx_inv = 1.0 / (double)nx;
    for (int i=0; i<=nx; ++i) {
        double x = x0 + (x1-x0) * (double)i * nx_inv;
        spec.xb.push_back(x);
    }

    // Set up y coordinates
    int ny = (int)(.5 + (y1 - y0) / dy);    // Round to nearest integer
    double ny_inv = 1.0 / (double)ny;
    for (int i=0; i<=ny; ++i) {
        double y = y0 + (y1-y0) * (double)i * ny_inv;
        spec.yb.push_back(y);
    }

    return spec;
}

void Grid make_grid(
    std::string const &name,
    GridSpec_XY const &spec,
    std::function<bool(Cell const &)> euclidian_clip = &EuclidianClip::keep_all)
{
    auto &xb(spec.xb);
    auto &yb(spec.yb);

    Indexing indexing({"x", "y"}, {0,0}, {spec.nx(), spec.ny()}, spec.indices);

    // Set up the main grid
    GridMap<Vertex> vertices(xb.size() * yb.size());
    GridMap<Cell> cells(spec.nx() * spec.ny());

    VertexCache vcache(&vertices);
    for (int iy = 0; iy < yb.size()-1; ++iy) {      // j
        double y0 = yb[iy];
        double y1 = yb[iy+1];

        for (int ix = 0; ix < xb.size()-1; ++ix) {      // i
            long index = indexing.tuple_to_index<int,2>({ix, iy});

            double x0 = xb[ix];
            double x1 = xb[ix+1];

            Cell cell;
            vcache.add_vertex(cell, x0, y0);
            vcache.add_vertex(cell, x1, y0);
            vcache.add_vertex(cell, x1, y1);
            vcache.add_vertex(cell, x0, y1);

            // Don't include things outside our clipping region
            if (!euclidian_clip(cell)) continue;

            cell.index = index;
            cell.i = ix;
            cell.j = iy;
            cell.native_area = cell.proj_area(nullptr);

            cells.add(std::move(cell));
        }
    }

    return Grid(name, GridType::XY,
        GridCoordinates::XY, sproj,
        Gridparameterization::L0,
        std::move(indexing),
        std::unique_ptr<GridSpec>(new GridSpec_XY(spec)),
        std::move(vertices), std::move(cells));
}

// ---------------------------------------------------------



}
