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

#include <icebin/error.hpp>
#include <icebin/gridgen/gridutil.hpp>
#include <icebin/gridgen/clippers.hpp>
#include <icebin/gridgen/GridSpec_XY.hpp>

using namespace std::placeholders;  // for _1, _2, _3...
using namespace ibmisc;
using namespace icebin;
namespace po = boost::program_options;


static const double km = 1000.0;

BOOST_ENUM_VALUES(IceModel, int,
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
        ("icemodel", po::value<std::string>(), "Ice model to use (pism or searise)");

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

    int grid_size = 25;
    if (vm.count("grid")) grid_size = vm["grid"].as<int>();

    IceModel icemodel = IceModel::pism;
    if (vm.count("icemodel")) {
        std::string sicemodel = vm["icemodel"].as<std::string>();
        icemodel = ibmisc::parse_enum<IceModel>(sicemodel);
    }

    // ------------------------------------------------------
    double dsize = (double)grid_size;

    printf("------------- Set up the local ice grid\n");

    GridSpec_XY spec;
    char xname[100];
    snprintf(xname, sizeof(xname), "mar_g%d_%s", grid_size, icemodel.str());
    std::string name(xname);

    spec.name = name;
    spec.euclidian_clip = &EuclidianClip::keep_all;

    if (zone == Zone::greenland) {
        // The true exact MAR grid
        //spec.sproj = "+proj=stere +lon_0=-40 +lat_0=70.5 +lat_ts=0 +ellps=WGS84";
	//set_xy_boundaries(spec,
        //    (- 800.0 - .5*dsize)*km, (- 800.0 + 300.0*5 + .5*dsize)*km,   dsize*km,
        //    (-3400.0 - .5*dsize)*km, (-3400.0 + 560.0*5 + .5*dsize)*km,   dsize*km);
        spec.sproj = "+proj=stere +lon_0=-40.0 +lat_0=70.5 +lat_ts=0.0 +k=1.0 +a=6371229 +b=6371129 +no_defs";
        set_xy_boundaries(spec,
            (- 775.0 - .5*dsize)*km, (700 + .5*dsize)*km,   dsize*km,
            (-1200.0 - .5*dsize)*km, (1525 + .5*dsize)*km,   dsize*km);
	//For MAR 25 km version only.  Different resolution runs may have different domain sizes.
	//
    } else {    // Antarctica
        // Zone must be Greenland, for now...
        //spec.sproj = "+proj=stere +lon_0=0 +lat_0=-90 +lat_ts=71.0 +ellps=WGS84";

        // The true exact MAR grid
        //set_xy_boundaries(spec,
        //    (-2800.0 - .5*dsize)*km, (-2800.0 + 1200*5.0 + .5*dsize)*km,   dsize*km,
        //   (-2800.0 - .5*dsize)*km, (-2800.0 + 1200*5.0 + .5*dsize)*km,   dsize*km);
    }


    // ------------ Make the grid from the spec
    Grid_XY grid;
    if (icemodel == IceModel::pism) {
        spec.indexing = Indexing({"x", "y"}, {0,0}, {spec.nx(), spec.ny()}, {0,1});   // row major: x has largest stride
    } else {    // Native MAR (and ModelE)
        spec.indexing = Indexing({"x", "y"}, {0,0}, {spec.nx(), spec.ny()}, {1,0});   // column major: y has largest stride
    }
    spec.make_grid(grid);

    // ------------- Write it out to NetCDF
    ibmisc::NcIO ncio(spec.name + ".nc", netCDF::NcFile::replace);
    grid.ncio(ncio, "grid");
    ncio.close();
}
