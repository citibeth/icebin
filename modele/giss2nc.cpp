// Convert GISS-format binary Fortran files to NetCDF

#include <iostream>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <ibmisc/fortranio.hpp>
#include <ibmisc/memory.hpp>
#include <ibmisc/blitz.hpp>
#include <icebin/error.hpp>
#include <icebin/modele/z1qx1n_bs1.hpp>

using namespace std;
using namespace ibmisc;
using namespace icebin;
using namespace icebin::modele;
namespace po = boost::program_options;




std::vector<fortran::Shape<2>> stdshapes {
    fortran::Shape<2>({"im2", "jm2"}, {IM2, JM2}),
    fortran::Shape<2>({"ims", "jms"}, {IMS, JMS}),
    fortran::Shape<2>({"imh", "jmh"}, {IMH, JMH}),
    fortran::Shape<2>({"im1", "jm1"}, {IM1, JM1}),
    fortran::Shape<2>({"im", "jm"}, {IM, JM})
};

// ----------------------------------------------------------------------
static const boost::regex
    titleiRE("(.*?):(.*?)(\\((.*?)\\))?\\s\\s(.*)");

struct Info {
    std::string name;
    std::string description;
    std::string units;
    std::string source;

    Info(
        std::string const &_name,
        std::string const &_description,
        std::string const &_units,
        std::string const &_source)
    : name(_name), description(_description), units(_units), source(_source) {}
};

std::string const get_match(boost::cmatch const &matches, int i)
{
    std::string str(matches[i].str());
    boost::trim(str);
    return str;
}

Info parse_titlei(std::string const &str)
{
    boost::cmatch matches;
    if (boost::regex_match(str.c_str(), matches, titleiRE)) {
        return Info(
            get_match(matches, 1),
            get_match(matches, 2),
            get_match(matches, 4),
            get_match(matches, 5));
    } else {
        (*icebin_error)(-1,
            "Cannot parse titlei='%s'", str.c_str());
    }
}
// ----------------------------------------------------------------------

int main(int argc, char **argv)
{
    // Parse comamnd line arguments
    std::string ifname, ofname;
    try {
	    po::options_description desc("Allowed options");
	    desc.add_options()
	        ("help,h", "produce help message")
	        ("input-file", po::value<string>(), "input file")
	        ("output-file", po::value<string>(), "output file")
	    ;

	    po::positional_options_description positional;
	    positional.add("input-file", 1);
	    positional.add("output-file", 1);

	    if (argc == 1) {
	        cout << desc << endl;
	        return 1;
	    }

	    po::variables_map vm;
	    po::command_line_parser parser(argc, argv);
	    po::store(po::command_line_parser(argc, argv).
	          options(desc).positional(positional).run(), vm);
	    po::notify(vm);
	    
	    if (vm.count("help")) {
	        cerr << desc << endl;
	        return 1;
	    }

        ifname = vm["input-file"].as<std::string>();
        ofname = vm["output-file"].as<std::string>();

    } catch(std::exception &exp) {
        cout << "Error parsing arumgnets:" << endl << exp.what() << endl;
        return 1;
    }

    printf("ARGS: %s %s\n", ifname.c_str(), ofname.c_str());


    std::array<char, 80> titlei;
    fortran::Shape<2> const *data_shape;

    {fortran::UnformattedInput fin(ifname, Endian::LITTLE);
    TmpAlloc tmp;
    NcIO ncio(ofname, 'w');

        // Read and then write a single array
        for (;;) {
	        auto &data(tmp.make<blitz::Array<float,2>>());

            // Read from Fortran binary file
	        fortran::read(fin) >> titlei
	            >> fortran::star(data, data_shape, stdshapes)
	            >> fortran::endr;

            // EOF gets set if we tried to read off the end of the file
            if (fin.eof()) break;

            // Parse titlei
            // https://panthema.net/2007/0314-BoostRegex
            Info info(parse_titlei(fortran::trim(titlei)));

            // Write to NetCDF
	        auto ncvar(ncio_blitz(ncio, data, false, info.name,
	            get_or_add_dims(ncio,
                    to_vector(data_shape->sshape),
                    to_vector_cast<int,long,2>(data_shape->shape))));
            get_or_put_att(ncvar, ncio.rw, "description", info.description);
            get_or_put_att(ncvar, ncio.rw, "units", info.units);
            get_or_put_att(ncvar, ncio.rw, "source", info.source);
        }
    }
}
