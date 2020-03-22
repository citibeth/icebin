#include <iostream>
#include <tclap/CmdLine.h>
#include <icebin/modele/topo_base.hpp>
#include <ibmisc/memory.hpp>
#include <everytrace.h>

using namespace ibmisc;
using namespace icebin::modele;

struct ParseArgs {
    std::string ofname;
    std::string et1mfile;

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

        TCLAP::UnlabeledValueArg<std::string> et1mfile_a(
            "et1mfile", "IN: Name of ETOPO1 file", true, "", "ETOPO1 filename", cmd);

        TCLAP::UnlabeledValueArg<std::string> ofname_a(
            "ofname", "OUT: Name of output file", true, "", "output filename", cmd);

        // Parse the argv array.
        cmd.parse( argc, argv );

        // Get the value parsed by each arg.
        et1mfile = et1mfile_a.getValue();
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

    MakeTopoO mto(EnvSearchPath("MODELE_FILE_PATH"), {
        "FGICE1m", args.et1mfile, "FGICE1m",
        "ZICETOP1m", args.et1mfile, "ZICETOP1m",
        "ZSOLG1m", args.et1mfile, "ZSOLG1m",
        "FOCEAN1m", args.et1mfile, "FOCEAN1m",
        "FLAKES", "Z10MX10M.nc", "FLAKES"
    });


    // Print sanity check errors to STDERR
    for (std::string const &err : mto.errors) fprintf(stderr, "ERROR: %s\n", err.c_str());


    printf("============ Writing Output\n");
    std::vector<double> lonc(mto.hspec.lonc());
    std::vector<double> latc(mto.hspec.latc());
    {NcIO ncio(args.ofname, 'w');
        ncio_vector(ncio, lonc, false, "lon", "double",
            get_or_add_dims(ncio, {"im"}, {(long)lonc.size()}));
        ncio_vector(ncio, latc, false, "lat", "double",
            get_or_add_dims(ncio, {"jm"}, {(long)latc.size()}));

        mto.hspec.ncio(ncio, "hspec");
        mto.bundle.ncio(ncio, {}, "", "double");

        ncio.close();
    }
    printf("Done Writing Output\n");

    if (mto.errors.size() > 0) return -1;
    return 0;
}
