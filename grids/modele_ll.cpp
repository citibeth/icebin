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

#include <ibmisc/indexing.hpp>

#include <icebin/error.hpp>
#include <icebin/gridgen/gridutil.hpp>
#include <icebin/gridgen/GridSpec_LonLat.hpp>
#include <icebin/modele/GridSpec_Hntr.hpp>
#include <icebin/modele/clippers.hpp>
#include <icebin/modele/grids.hpp>

using namespace std::placeholders;  // for _1, _2, _3...
using namespace ibmisc;
using namespace icebin;
using namespace icebin::modele;
namespace po = boost::program_options;

const double km = 1000.0;
const double EQ_RAD = 6.371e6; /// Radius of the Earth (same as in ModelE)


// http://stackoverflow.com/questions/313970/how-to-convert-stdstring-to-lower-case
void tolower(std::string &data)
{
    for (auto ii(data.begin()); ii != data.end(); ++ii)
        *ii = std::tolower(*ii);
}

std::map<std::string, HntrGrid const *> grids_by_name {
    {"2mx2m", &g2mx2m},
    {"10mx10m", &g10mx10m},
    {"hxh", &ghxh},
    {"1x1", &g1x1},
    {"1qx1", &g1qx1},
    {"2hx2", &g2hx2},
    {"2x2_5", &g2hx2},    // Legacy compatibility
    {"5x4", &g5x4},
};

int main(int argc, char **argv)
{
    ibmisc::netcdf_debug = true;

    // -------------------------------------------------------------
    // Parse Command Line Args

    // http://www.boost.org/doc/libs/1_60_0/doc/html/program_options/tutorial.html
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("zone", po::value<std::string>(), "Greenland, Antarctica or Both")
        ("grid", po::value<std::string>(), "2mx2m, 10mx10m, hxh, 1x1, q1x1, 2hx2, 5x4")
        ("pole-caps", po::value<bool>()->default_value(true), "Combined cap grid cell at poles?")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (argc == 0 || vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    int zone;
    std::string szone = "ga";
    if (vm.count("zone")) szone = vm["zone"].as<std::string>();
    tolower(szone);
    if (szone == "greenland" || szone == "g") {
        zone = ice_sheet::GREENLAND;
        szone = "g";
    } else if (szone == "antarctica" || szone == "a") {
        zone = ice_sheet::ANTARCTICA;
        szone = "a";
    } else if (szone == "both" || szone == "ga") {
        zone = ice_sheet::GREENLAND | ice_sheet::ANTARCTICA;
        szone = "ga";
    } else {
        (*icebin_error)(-1, "Unknown zone '%s'\n", szone.c_str());
    }

    std::string sgrid = "2hx2";
    if (vm.count("grid")) sgrid = vm["grid"].as<std::string>();
    HntrGrid const *hntr_grid = grids_by_name.at(sgrid);

    // -------------------------------------------------------------
    GridSpec_Hntr spec(*hntr_grid);
    spec.name = "modele_ll_" + szone + sgrid;
    spec.pole_caps = vm["pole-caps"].as<bool>();

    spec.spherical_clip = std::bind(&ice_sheet::clip, zone, _1, _2, _3, _4, _5);
    spec.points_in_side = (hntr_grid->im > IM1 ? 1 : 2);    // Use 2 for comparison with past experiments
    spec.eq_rad = EQ_RAD;
    

    // ------------ Make the grid from the spec
    Grid_LonLat grid;
    spec.make_grid(grid);

    // ------------- Write it out to NetCDF
    ibmisc::NcIO ncio(spec.name + ".nc", netCDF::NcFile::replace);
    grid.ncio(ncio, "grid");
    ncio.close();
}
