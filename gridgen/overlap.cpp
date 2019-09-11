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
#include <tclap/CmdLine.h>

#include <ibmisc/netcdf.hpp>
#include <ibmisc/stdio.hpp>

#include <icebin/Grid.hpp>
#include <icebin/gridgen/GridGen_Exchange.hpp>
#include <icebin/gridgen/gridutil.hpp>

static const double km = 1000.0;

using namespace icebin;
using namespace ibmisc;
using namespace netCDF;

struct ParseArgs {
    std::string fnameA;
    std::string fnameI;
    std::string fname_exgrid;    // OUT: Name of overlap file to write

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

        TCLAP::UnlabeledValueArg<std::string> fnameA_a("gridA",
            "Name of file containing Atmosphere grid",
            true, "", "atmosphere grid file", cmd);

        TCLAP::UnlabeledValueArg<std::string> fnameI_a("gridI",
            "Name of file containing Ice grid",
            true, "", "ice grid file", cmd);


        TCLAP::ValueArg<std::string> fname_exgrid_a("o", "out",
            "Name of IceBin overlap file to write",
            false, "", "overlap grid file", cmd);


        // Parse the argv array.
        cmd.parse( argc, argv );

        fnameA = fnameA_a.getValue();
        fnameI = fnameI_a.getValue();
        fname_exgrid = fname_exgrid_a.getValue();
    } catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        exit(1);
    }
}


int main(int argc, char **argv)
{
    ParseArgs args(argc, argv);

    printf("------------- Read gridA (GCM Grid): %s\n", args.fnameA.c_str());
    NcIO ncio1(args.fnameA, 'r');
    Grid gridA;
    gridA.ncio(ncio1, "grid");
    ncio1.close();
    printf("Done reading gridA\n");

    printf("------------- Read gridI (Ice Grid): %s\n", args.fnameI.c_str());
    NcIO ncio2(args.fnameI, 'r');
    Grid gridI;
    gridI.ncio(ncio2, "grid");
    ncio2.close();
    printf("Done reading gridI\n");

    printf("--------------- Overlapping\n");
    Grid exgrid(make_exchange_grid(&gridA, &gridI));
    sort_renumber_vertices(exgrid);

    printf("--------------- Writing out\n");
    std::string fname(args.fname_exgrid);
    if (fname == "")
        fname = strprintf("%s.nc", exgrid.name.c_str());    // Using operator+() or append() doesn't work here with GCC 4.9.3

    printf("overlap writing to %s", fname.c_str());
    ibmisc::NcIO ncio(fname, 'w');
    gridA.ncio(ncio, "gridA");
    gridI.ncio(ncio, "gridI");
    exgrid.ncio(ncio, "exgrid");
    ncio.close();
}
