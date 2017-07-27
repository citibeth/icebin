#include <iostream>
#include <tclap/CmdLine.h>
#include <icebin/modele/z1qx1n_bs1.hpp>
#include <ibmisc/memory.hpp>
#include <everytrace.h>

using namespace ibmisc;
using namespace icebin::modele;

struct ParseArgs {
    std::string ofname;

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


        // Parse the argv array.
        cmd.parse( argc, argv );

        // Get the value parsed by each arg.
        ofname = ofname_a.getValue();

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

    // Read the input files
    printf("============ Reading Input Files\n");
    TopoInputs topo_inputs(true);
    read_raw(topo_inputs, NULL, EnvSearchPath("MODELE_FILE_PATH"));

    // Do the calculation
    printf("============ Calculating TOPO\n");
    TopoOutputs topo_outputs(true);
    z1qx1n_bs1(topo_inputs, topo_outputs);

    printf("============ Writing Output\n");
    {NcIO ncio(args.ofname, 'w');

        topo_outputs.bundle.ncio(ncio, {}, false, "", "double");

        ncio.close();
    }
    printf("Done Writing Output\n");

}
