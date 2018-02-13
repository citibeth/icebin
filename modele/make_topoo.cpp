#include <iostream>
#include <tclap/CmdLine.h>
#include <icebin/modele/make_topoo.hpp>
#include <ibmisc/memory.hpp>
#include <everytrace.h>

using namespace ibmisc;
using namespace icebin::modele;

struct ParseArgs {
    std::string ofname;
    bool greenland;

    ParseArgs(int argc, char **argv);
};

ParseArgs::ParseArgs(int argc, char **argv)
{
    // Wrap everything in a try block.  Do this every time, 
    // because exceptions will be thrown for problems.
    try {  
        TmpAlloc tmp;

        // Define the command line object, and insert a message
        // that describes the program. The "Command description message" 
        // is printed last in the help text. The second argument is the 
        // delimiter (usually space) and the last one is the version number. 
        // The CmdLine object parses the argv array based on the Arg objects
        // that it contains. 
        TCLAP::CmdLine cmd("Command description message", ' ', "<no-version>");

        TCLAP::UnlabeledValueArg<std::string> ofname_a(
            "ofname", "Name of output file", true, "", "output filename", cmd);
        TCLAP::SwitchArg greenland_a("g", "greenland", "Include Greenland?", cmd, false);

        // Parse the argv array.
        cmd.parse( argc, argv );

        // Get the value parsed by each arg.
        ofname = ofname_a.getValue();
        greenland = greenland_a.getValue();

    } catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        exit(1);
    }
}



int main(int argc, char** argv)
{
    using namespace TCLAP;

    everytrace_init();
    ParseArgs args(argc, argv);

    std::string et1mfile("etopo1_ice_g1m.nc");
    make_topoO(EnvSearchPath("MODELE_FILE_PATH"), {
        "FGICE1m", et1mfile, "FGICE1m",
        "ZICETOP1m", et1mfile, "ZICETOP1m",
        "ZSOLG1m", et1mfile, "ZSOLG1m",
        "FOCEAN1m", et1mfile, "FOCEAN1m",
        "FLAKES", "Z10MX10M.nc", "FLAKES"
    }));


    printf("============ Writing Output\n");
    {NcIO ncio(args.ofname, 'w');

        bundle.ncio(ncio, {}, "", "double");

        ncio.close();
    }
    printf("Done Writing Output\n");

}
