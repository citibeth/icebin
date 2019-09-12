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

// Generates PISM grid as of year 2019

#include <string>
#include <ctype.h>
#include <algorithm>
#include <functional>
#include <iostream>

#include <boost/filesystem.hpp>
#include <tclap/CmdLine.h>

#include <ibmisc/enum.hpp>
#include <ibmisc/stdio.hpp>

#include <icebin/error.hpp>
#include <icebin/gridgen/gridutil.hpp>
#include <icebin/gridgen/clippers.hpp>
#include <icebin/gridgen/GridGen_XY.hpp>
#include <icebin/gridgen/GridGen_LonLat.hpp>

using namespace ibmisc;
using namespace icebin;
using namespace netCDF;

struct ParseArgs {
    std::string spec;    // INPUT file
    std::string grid;    // OUTPUT file

    ParseArgs(int argc, char **argv);
};

ParseArgs::ParseArgs(int argc, char **argv)
{
    // Wrap everything in a try block.  Do this every time, 
    // because exceptions will be thrown for problems.
    try {  
        // Define the command line object, and insert a message
        // that describes the program. The "Command description message" 
        // is printed last in the help text. The second argument is the 
        // delimiter (usually space) and the last one is the version number. 
        // The CmdLine object parses the argv array based on the Arg objects
        // that it contains. 
        TCLAP::CmdLine cmd("Command description message", ' ', "<no-version>");

        TCLAP::UnlabeledValueArg<std::string> spec_a("input",
            "Name of PISM state file, after PISM spinup",
            true, "", "pism state file", cmd);

        TCLAP::ValueArg<std::string> grid_a("o", "output",
            "Name of IceBin grid file to write",
            false, "", "icebin grid file", cmd);

        // Parse the argv array.
        cmd.parse( argc, argv );

        spec = spec_a.getValue();
        grid = grid_a.getValue();
    } catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        exit(1);
    }
}

int main(int argc, char **argv)
{
    // Parse Command Line Args
    ParseArgs args(argc,argv);

    std::unique_ptr<GridSpec> spec;
    std::string name;
    {NcIO ncio(args.spec, 'r');
        NcVar info_v = get_or_add_var(ncio, "grid.info", "int", {});
        get_or_put_att(info_v, ncio.rw, "name", name);
        ncio_grid_spec(ncio, spec, "grid");
    }
 
    // ------------ Make the grid from the spec
    std::unique_ptr<Grid> grid;
    switch(spec->type.index()) {
        case GridType::XY :
            grid.reset(new Grid(make_grid(name, *dynamic_cast<GridSpec_XY *>(&*spec), &EuclidianClip::keep_all)));
        break;
        case GridType::LONLAT :
            grid.reset(new Grid(make_grid(name, *dynamic_cast<GridSpec_LonLat *>(&*spec), &SphericalClip::keep_all)));
        break;
        default:
            fprintf(stderr, "Do not know how to generate grid type of type %s", spec->type.str());
            return -1;
    }

    // ------------- Write it out to NetCDF
    std::string ofname(args.grid);
    if (ofname == "") ofname = strprintf("%s.nc", name.c_str());    // Using operator+() or append() doesn't work here with GCC 4.9.3

    ibmisc::NcIO ncio(ofname, netCDF::NcFile::replace);
    grid->ncio(ncio, "grid");
    ncio.close();
}
