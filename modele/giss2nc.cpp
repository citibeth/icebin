// Convert GISS-format binary Fortran files to NetCDF

#include <functional>
#include <iostream>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <ibmisc/fortranio.hpp>
#include <ibmisc/memory.hpp>
#include <ibmisc/blitz.hpp>
#include <ibmisc/endian.hpp>
#include <ibmisc/string.hpp>

#include <icebin/error.hpp>
#include <icebin/modele/grids.hpp>

using namespace std;
using namespace ibmisc;
using namespace icebin;
using namespace icebin::modele;
namespace po = boost::program_options;


std::vector<fortran::Shape<2>> stdshapes {
    fortran::Shape<2>({"im1m", "jm1m"}, {IM1m, JM1m}),
    fortran::Shape<2>({"im1m", "jm1m_gridreg"}, {IM1m, JM1m_gridreg}),
    fortran::Shape<2>({"im2m", "jm2m"}, {IM2m, JM2m}),
    fortran::Shape<2>({"ims", "jms"}, {IMS, JMS}),
    fortran::Shape<2>({"imh", "jmh"}, {IMH, JMH}),
    fortran::Shape<2>({"im1", "jm1"}, {IM1, JM1}),
    fortran::Shape<2>({"im", "jm"}, {IM, JM}),
    fortran::Shape<2>({"im5m", "jm5m"}, {IM5m, JM5m})
};

// ----------------------------------------------------------------------
static const boost::regex
    titleiRE("(.*?):(.*?)(\\((.*?)\\))?(\\s\\s(.*))?");

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

Info parse_titlei(std::string const &str, std::string const &use_name)
{
    boost::cmatch matches;
    if (boost::regex_match(str.c_str(), matches, titleiRE)) {
        if (use_name == "") {
            return Info(
                get_match(matches, 1),
                get_match(matches, 2),
                get_match(matches, 4),
                get_match(matches, 6));
        } else {
            return Info(
                use_name,
                get_match(matches, 1) + " " + get_match(matches, 2),
                get_match(matches, 4),
                get_match(matches, 6));
        }
    } else {
        return Info(use_name, str, "", "");
    }
}
// ----------------------------------------------------------------------
template<class IType, class OType>
void giss2nc(
    std::string const &ifname,
    std::string const &ofname,
    Endian endian,
    std::vector<std::string> const &names)
{
    std::array<char, 80> titlei;
    fortran::Shape<2> const *data_shape;

    {fortran::UnformattedInput fin(ifname, endian);
    TmpAlloc tmp;
    NcIO ncio(ofname, 'w');

        // Read and then write a single array
        for (size_t i=0;; ++i) {

            if (names.size() > 0 && i >= names.size()) break;

	        auto &data(tmp.make<blitz::Array<IType,2>>());

            // Read from Fortran binary file
            // Allocates as blitz::fortranArray
	        fortran::read(fin) >> titlei
	            >> fortran::allocate(data, data_shape, stdshapes)
	            >> fortran::endr;

            // EOF gets set if we tried to read off the end of the file
            if (fin.eof()) break;

            // Parse titlei
            // https://panthema.net/2007/0314-BoostRegex
            std::string const &use_name = (i < names.size() ? names[i] : "");
            std::string titlei_cxx(fortran::trim(titlei));
            Info info(parse_titlei(titlei_cxx, use_name));
            if (info.name == "") {
                fprintf(stderr,
                    "Cannot parse titlei='%s'; try using names parameter\n", titlei_cxx.c_str());
                info.name = string_printf("var%02d", i);
            }

            printf("Read variable named %s: description=\"%s\"\n", info.name.c_str(), info.description.c_str());
            if (info.name != "_") {

                // Write to NetCDF
	            auto ncvar(ncio_blitz(
                    ncio, data, info.name, get_nc_type<OType>(),
    	            get_or_add_dims(ncio,
                        to_vector(data_shape->sshape),
                        to_vector_cast<int,long,2>(data_shape->shape)),
                    DimOrderMatch::MEMORY, false));
                get_or_put_att(ncvar, ncio.rw, "description", info.description);
                get_or_put_att(ncvar, ncio.rw, "units", info.units);
                get_or_put_att(ncvar, ncio.rw, "source", info.source);
            }
        }
    }
}
// ----------------------------------------------------------------------
// Build multiple instantiations of giss2nc(), let user choose at runtime.
typedef std::map<std::string,
    std::function<void(
        std::string const &ifname,
        std::string const &ofname,
        Endian endian,
        std::vector<std::string> const &names)>>
Giss2NcMap;

Giss2NcMap make_giss2nc_map()
{
    Giss2NcMap ret;
    ret.insert(std::make_pair("float-double", &giss2nc<float,double>));
    ret.insert(std::make_pair("int16", &giss2nc<int16_t,int16_t>));
    return ret;
}
Giss2NcMap giss2nc_map(make_giss2nc_map());

// ----------------------------------------------------------------------
int main(int argc, char **argv)
{
    // Parse comamnd line arguments
    std::string ifname, ofname;
    Endian endian;
    std::vector<std::string> names;
    std::string stype;

    try {
	    po::options_description desc("Allowed options");
	    desc.add_options()
	        ("help,h", "produce help message")
	        ("input-file", po::value<string>(), "input file")
	        ("output-file", po::value<string>(), "output file")
            ("endian", po::value<Endian>()->default_value(Endian::LITTLE), "[big | little]")
            ("type", po::value<string>()->default_value("float-double"),
                "Data type of (a) Fortran input, (b) NetCDF output [float-double | int16]")
            ("names", po::value<string>()->default_value(""), "(OPTIONAL) comma-separated variable names")
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
        endian = vm["endian"].as<Endian>();
        stype = vm["type"].as<std::string>();

        std::string snames = vm["names"].as<std::string>();
        if (snames != "")
            boost::algorithm::split(names, snames, boost::is_any_of(","));

    } catch(std::exception &exp) {
        cout << "Error parsing arumgnets:" << endl << exp.what() << endl;
        return 1;
    }

    printf("ARGS: %s %s\n", ifname.c_str(), ofname.c_str());

    auto _giss2nc(giss2nc_map.at(stype));
    _giss2nc(ifname, ofname, endian, names);
}

