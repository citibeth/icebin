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
#include <ctype.h>
#include <algorithm>
#include <functional>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <ibmisc/enum.hpp>
#include <ibmisc/stdio.hpp>

#include <icebin/error.hpp>
#include <icebin/gridgen/gridutil.hpp>
#include <icebin/gridgen/clippers.hpp>
#include <icebin/gridgen/GridGen_XY.hpp>

using namespace std::placeholders;  // for _1, _2, _3...
using namespace ibmisc;
using namespace icebin;
namespace po = boost::program_options;


static const double km = 1000.0;

BOOST_ENUM_VALUES(IndexOrder, int,
    (pism) (0)
    (searise) (1)
)

BOOST_ENUM_VALUES(Zone, int,
    (greenland) (0)
    (antarctica) (1)
)

int main(int argc, char **argv)
{
    // -------------------------------------------------------------
    // Parse Command Line Args

    // http://www.boost.org/doc/libs/1_60_0/doc/html/program_options/tutorial.html
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("zone", po::value<std::string>(), "Greenland or Antarctica")
        ("grid", po::value<int>(), "Cell size [km]")
        ("out,o", po::value<std::string>(), "Greenland or Antarctica")
        ("index_order", po::value<std::string>(), "Ice model to use (pism or searise)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (argc == 0 || vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    Zone zone = Zone::greenland;
    if (vm.count("zone")) {
        std::string szone = vm["zone"].as<std::string>();
        zone = ibmisc::parse_enum<Zone>(szone);
    }

    int grid_size = 20;
    if (vm.count("grid")) grid_size = vm["grid"].as<int>();

    std::string ofname;
    if (vm.count("out")) ofname = vm["out"].as<std::string>();

    IndexOrder index_order = IndexOrder::pism;
    if (vm.count("index_order")) {
        std::string sindex_order = vm["index_order"].as<std::string>();
        index_order = ibmisc::parse_enum<IndexOrder>(sindex_order);
    }


    // ------------------------------------------------------
    double dsize = (double)grid_size;

    printf("------------- Set up the local ice grid\n");

    // Get the spec
    // NOTE: PISM option is only for the OLD PISM.
    //       More recently, PISM switched their indices to be the same
    //       as SeaRISE and ModelE
    auto indices(index_order == IndexOrder::pism
        ? std::vector<int>{0,1}
        : std::vector<int>{1,0});

    GridSpec_XY spec;
    if (zone == Zone::greenland) {
        spec = GridSpec_XY::make_with_boundaries(
            "+proj=stere +lon_0=-39 +lat_0=90 +lat_ts=71.0 +ellps=WGS84",
            std::move(indices),
            (- 800.0 - .5*dsize)*km, (- 800.0 + 300.0*5 + .5*dsize)*km,   dsize*km,
            (-3400.0 - .5*dsize)*km, (-3400.0 + 560.0*5 + .5*dsize)*km,   dsize*km);
    } else {    // Antarctica
        spec = GridSpec_XY::make_with_boundaries(
            "+proj=stere +lon_0=0 +lat_0=-90 +lat_ts=71.0 +ellps=WGS84",
            std::move(indices),
            (-2800.0 - .5*dsize)*km, (-2800.0 + 1200*5.0 + .5*dsize)*km,   dsize*km,
            (-2800.0 - .5*dsize)*km, (-2800.0 + 1200*5.0 + .5*dsize)*km,   dsize*km);
    }

    // ------------ Make the grid from the spec
    std::string name(ibmisc::strprintf("sr_g%d_%s", grid_size, index_order.str()));
    Grid grid(make_grid(name, spec, &EuclidianClip::keep_all));

    // ------------- Write it out to NetCDF
    if (ofname == "") ofname = strprintf("%s.nc", name.c_str());    // Using operator+() or append() doesn't work here with GCC 4.9.3

    ibmisc::NcIO ncio(ofname, netCDF::NcFile::replace);
    grid.ncio(ncio, "grid");
    ncio.close();
}
