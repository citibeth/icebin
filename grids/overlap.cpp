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

#include <string>
#include <boost/filesystem.hpp>

#include <ibmisc/netcdf.hpp>

#include <icebin/Grid.hpp>
#include <icebin/gridgen/GridSpec_Exchange.hpp>
#include <icebin/gridgen/gridutil.hpp>

static const double km = 1000.0;

using namespace icebin;
using namespace ibmisc;
using namespace netCDF;

int main(int argc, char **argv)
{
    std::string fnameA(argv[1]);
    std::string fnameI(argv[2]);

    printf("------------- Read gridA (GCM Grid): %s\n", fnameA.c_str());
    NcIO ncio1(fnameA, 'r');
    Grid gridA;
    gridA.ncio(ncio1, "grid");
    ncio1.close();
    printf("Done reading gridA\n");

    printf("------------- Read gridI (Ice Grid): %s\n", fnameI.c_str());
    NcIO ncio2(fnameI, 'r');
    Grid gridI;
    gridI.ncio(ncio2, "grid");
    ncio2.close();
    printf("Done reading gridI\n");

    printf("--------------- Overlapping\n");
    GridSpec_Exchange spec;
    spec.gridA = &gridA;
    spec.gridI = &gridI;
    Grid exgrid;
    spec.make_grid(exgrid);
    sort_renumber_vertices(exgrid);

    printf("--------------- Writing out\n");
    std::string fname = exgrid.name + ".nc";

    ibmisc::NcIO ncio(exgrid.name + ".nc", 'w');
    gridA.ncio(ncio, "gridA");
    gridI.ncio(ncio, "gridI");
    exgrid.ncio(ncio, "exgrid");
    ncio.close();
}
