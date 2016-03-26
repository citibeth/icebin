/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <boost/bind.hpp>
#include <glint2/Grid_LonLat.hpp>
#include <glint2/clippers.hpp>
#include <boost/filesystem.hpp>
#include <glint2/modele/grids_ll.hpp>
#include <glint2/clippers.hpp>

static const double km = 1000.0;

using namespace glint2;

const double EQ_RAD = 6.371e6; /// Radius of the Earth (same as in ModelE)
//const double EQ_RAD = 6370997; /// Radius of the Earth (same as in proj.4, see src/pj_ellps.c)

int main(int argc, char **argv)
{
    boost::filesystem::path exe_path(argv[0]);
    printf("------------- Set up GCM Grid\n");

    Grid_LonLat grid;
    glint2::modele::set_lonlat_2x2_5(grid);
    grid.name = exe_path.stem().string();
//  grid.points_in_side = 4;
    grid.points_in_side = 2;    // Use 2 for comparison with past experiments
    grid.eq_rad = EQ_RAD;
    grid.realize(boost::bind(
        &SphericalClip::lonlat, -74., 59., -10., 87.5,
        _1, _2, _3, _4));

//  grid.realize(boost::bind(&SphericalClip::keep_all, _1, _2, _3, _4));

    // ------------- Write it out to NetCDF
    fflush(stdout);
    printf("// ------------- Write it out to NetCDF\n");
    grid.to_netcdf(grid.name + ".nc");
//  grid.to_netcdf("greenland_2x2_5.nc");
}
