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
#include <icebin/gridgen/clippers.hpp>
#include <icebin/gridgen/GridSpec_LonLat.hpp>

using namespace std::placeholders;  // for _1, _2, _3...
using namespace ibmisc;
using namespace icebin;
namespace po = boost::program_options;

const double km = 1000.0;
const double EQ_RAD = 6.371e6; /// Radius of the Earth (same as in ModelE)

const std::vector<double> lonb_4x5 = {-180,-175,-170,-165,-160,-155,-150,-145,-140,-135,-130,-125,-120,-115,-110,-105,-100,-95,-90,-85,-80,-75,-70,-65,-60,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180};

const std::vector<double> latb_4x5 = {-88,-84,-80,-76,-72,-68,-64,-60,-56,-52,-48,-44,-40,-36,-32,-28,-24,-20,-16,-12,-8,-4,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88};


// Bits...
const int GREENLAND = 1;
const int ANTARCTICA = 2;


/** Just initializes lonb and latb, doesn't do full set-up */
void set_lonlat_4x5(GridSpec_LonLat &spec)
{
    spec.lonb = lonb_4x5;
    spec.latb = latb_4x5;
    spec.south_pole = true;
    spec.north_pole = true;
}

void set_lonlat_2x2_5(GridSpec_LonLat &spec)
{
    // Create the 2x2.5 grid from the 4x5 grid.
    spec.lonb.clear();
    for (int i=0; i<lonb_4x5.size()-1; ++i) {
        spec.lonb.push_back(lonb_4x5[i]);
        spec.lonb.push_back(.5*(lonb_4x5[i] + lonb_4x5[i+1]));
    }
    spec.lonb.push_back(lonb_4x5.back());

    spec.latb.clear();
    for (int i=0; i<latb_4x5.size()-1; ++i) {
        spec.latb.push_back(latb_4x5[i]);
        spec.latb.push_back(.5*(latb_4x5[i] + latb_4x5[i+1]));
    }
    spec.latb.push_back(latb_4x5.back());

    spec.south_pole = true;
    spec.north_pole = true;
}

bool clip(int zone, double lon0, double lat0, double lon1, double lat1)
{
    // Is it in Greenland range?
    if (zone & GREENLAND)
        if (SphericalClip::lonlat(-74., 59., -10., 87.5,
            lon0, lat0, lon1, lat1)) return true;

    // Is it in Antarctica range?
    if (zone & ANTARCTICA)
        if (lat0 <= -60. || lat1 <= -60) return true;

    // Not in range of either ice sheet, discard
    return false;
}

// http://stackoverflow.com/questions/313970/how-to-convert-stdstring-to-lower-case
void tolower(std::string &data)
{
    for (auto ii(data.begin()); ii != data.end(); ++ii)
        *ii = std::tolower(*ii);
}

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
        ("grid", po::value<std::string>(), "4x5, 2x2.5");
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
        zone = GREENLAND;
        szone = "g";
    } else if (szone == "antarctica" || szone == "a") {
        zone = ANTARCTICA;
        szone = "a";
    } else if (szone == "both" || szone == "ga") {
        zone = GREENLAND | ANTARCTICA;
        szone = "ga";
    } else {
        (*icebin_error)(-1, "Unknown zone '%s'\n", szone.c_str());
    }

    std::string sgrid = "2x2_5";
    if (vm.count("grid")) sgrid = vm["grid"].as<std::string>();
    std::function<void (GridSpec_LonLat &)> spec_fn;
    if (sgrid == "2x2_5") {
         spec_fn = &set_lonlat_2x2_5;
    } else if (sgrid == "4x5") {
        spec_fn = &set_lonlat_4x5;
    } else {
        (*icebin_error)(-1, "Unknown grid '%s'\n", sgrid.c_str());
    }

    // -------------------------------------------------------------
    GridSpec_LonLat spec;
    spec.name = "modele_ll_" + szone + sgrid;

    spec.spherical_clip = std::bind(&clip, zone, _1, _2, _3, _4);
    spec_fn(spec);
    spec.points_in_side = 2;    // Use 2 for comparison with past experiments
    spec.eq_rad = EQ_RAD;
    

    // ------------ Make the grid from the spec
    Grid_LonLat grid;
    spec.indexing = Indexing<int,long>(
        {0,0}, {spec.nlon(), spec.nlat()}, {1,0});  // col major
    spec.make_grid(grid);

    // ------------- Write it out to NetCDF
    ibmisc::NcIO ncio(spec.name + ".nc", netCDF::NcFile::replace);
    grid.ncio(ncio, "grid");
    ncio.close();
}
