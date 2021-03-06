#include <cstdio>
#include <cstring>
#include <prettyprint.hpp>
#include <icebin/error.hpp>
#include <ibmisc/string.hpp>
#include <ibmisc/linear/linear.hpp>
#include <ibmisc/linear/compressed.hpp>
#include <icebin/modele/global_ec.hpp>
#include <tclap/CmdLine.h>

using namespace std;
using namespace ibmisc;
using namespace netCDF;
using namespace icebin;
using namespace icebin::modele;

static double const NaN = std::numeric_limits<double>::quiet_NaN();


std::array<std::string,2> parse_spec(std::string const &spec)
{
    auto dims(split<std::string>(spec, "v"));
    if (dims.size() != 2) (*icebin_error)(-1, "Error parsing matrix spec %s", spec.c_str());
    return std::array<std::string,2>{dims[0], dims[1]};
}

struct ParseArgs {
    bool scale;
    // Names of matrices to generate
    std::vector<std::array<std::string,2>> matrix_specs = vector<array<string,2>> {
        {"I","E"}, {"E","I"},
        {"I","A"}, {"A","I"},
        {"A","E"},
        {"I2", "A"}, {"I2", "E"}
    };

    std::vector<std::string> ifnames;    // Names of files to merge
    ParseArgs(int argc, char **argv);
};

ParseArgs::ParseArgs(int argc, char **argv)
{

    // Wrap everything in a try block.  Do this every time, 
    // because exceptions will be thrown for problems.
    try {  
        TCLAP::CmdLine cmd("Command description message", ' ', "<no-version>");

        TCLAP::ValueArg<std::string> matrix_names_a("n", "matrix-names",
            "Comma-separated names of matrices to generate, no spaces",
            false, "AvI,EvI,IvE,IvA,AvE", "matrix names", cmd);

        TCLAP::ValueArg<bool> scale_a("s", "scale",
            "Scale the matrix when joining?",
            false, false, "scale", cmd);

        TCLAP::UnlabeledMultiArg<std::string> ifnames_a("merge-files",
            "Files with sub-matrices to merge together",
            true, "filenames", cmd);

        cmd.parse(argc, argv);

        scale = scale_a.getValue();
        if (matrix_names_a.getValue() != "") {
            matrix_specs.clear();
            std::vector<std::string> matrix_names = split<std::string>(matrix_names_a.getValue(), ",");
            for (auto &name : matrix_names) matrix_specs.push_back(parse_spec(name));
        }

        ifnames = ifnames_a.getValue();
    } catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        exit(1);
    }
}



template<class T, int len>
blitz::TinyVector<T, len> vector_to_tiny(std::vector<T> const &vec)
{
    if (vec.size() != len) {
        (*ibmisc::ibmisc_error)(-1,
            "vector_to_tiny(): vector length %ld does not match declared length %d\n", vec.size(), len);
    }

    blitz::TinyVector<T, len> ret;
    for (int i=0; i < len; ++i) {
        ret[i] = vec[i];
    }
    return ret;
}

void combine_chunks(
    std::vector<std::string> const &ifnames,    // Names of input chunks
    std::string const &ofname,
    char ofmode,    // 'w' or 'a' for write mode of output file
    std::array<std::string,2> const &_sgrids,    // {"B", "A"} --> matrix BvA
    bool scale)
{


    printf("======== BEGIN combine_chunks(%sv%s)\n", _sgrids[0].c_str(), _sgrids[1].c_str());

    // Filled in below
    std::array<std::string,2> sgrids;
    std::string BvA;

    // Get total sizes
    int nnz = 0;        // = number non-zero]
    std::array<long,2> sparse_extents;
    std::array<std::vector<int>,2> shapes;
    std::vector<size_t> sizes;    // For printing

    global_ec::Metadata meta;
    for (size_t i=0; i<ifnames.size(); ++i) {
        std::string const &ifname(ifnames[i]);
        NcIO ncio(ifname, 'r');
        if (i == 0) {
            meta.ncio(ncio);

            // Names of the variables we will read/write
            // Use Ocean grid instead of Atmosphere grid if this is really an ocean file
            sgrids = _sgrids;
            std::replace(sgrids[0].begin(), sgrids[0].end(), 'A', 'O');
            std::replace(sgrids[1].begin(), sgrids[1].end(), 'A', 'O');
            BvA = sgrids[0] + "v" + sgrids[1];
        }

        netCDF::NcDim sz_nc = ncio.nc->getDim(BvA+".M.nnz");
        size_t sz = sz_nc.getSize();
        nnz += sz;
        sizes.push_back(sz);

        for (int i=0; i<2; ++i) {
            NcVar nc_dimA(ncio.nc->getVar("dim" + sgrids[i]));
            nc_dimA.getAtt("sparse_extent").getValues(&sparse_extents[i]);
            get_or_put_att(nc_dimA, 'r', "shape", "int", shapes[i]);
        }
    }
    cout << "Total non-zero elements in matrix = " << sizes << " = " << nnz << endl;
    cout << "sparse_extents = " << sparse_extents << endl;


    // -------- Get scale factors (sparse)

    long shape0;
    {NcIO ncio(ifnames[0], 'r');
        auto _dimB(ncio.nc->getVar("dim"+sgrids[0]));
        get_or_put_att(_dimB, 'r', "sparse_extent", "long", &shape0, 1);
    }
    blitz::Array<double,1> sM_s(shape0);
    sM_s = 0;

    // Start by accumulating weights, then invert to scale
    for (std::string const &ifname : ifnames) {
        printf("--- Reading %s\n", ifname.c_str());
        NcIO ncio(ifname, 'r');

        std::array<long,2> shape;
        auto _dimB(ncio.nc->getVar("dim"+sgrids[0]));
        get_or_put_att(_dimB, 'r', "sparse_extent", "long", &shape[0], 1);
        auto _dimA(ncio.nc->getVar("dim"+sgrids[1]));
        get_or_put_att(_dimA, 'r', "sparse_extent", "long", &shape[1], 1);

        auto dimB(nc_read_blitz<int,1>(ncio.nc, "dim"+sgrids[0]));
        auto wM_d(nc_read_blitz<double,1>(ncio.nc, BvA+".wM"));
        for (int i=0; i<wM_d.extent(0); ++i) {
            sM_s(dimB(i)) += wM_d(i);
        }
    }

    // Invert to get scale factors
    for (int i=0; i<sM_s.extent(0); ++i) sM_s(i) = 1. / sM_s(i);

    // -----------------------


    // Allocate
    linear::Weighted_Compressed ret;

    {auto wM(ret.weights[0].accum());
    auto M(ret.M.accum());
    auto Mw(ret.weights[1].accum());


        // Aggregate them together
        int iGlobal = 0;    // Count for debugging
        int iM=0;
        for (std::string const &ifname : ifnames) {
            printf("--- Reading %s\n", ifname.c_str());
            NcIO ncio(ifname, 'r');

            auto dimB(nc_read_blitz<int,1>(ncio.nc, "dim"+sgrids[0]));
            auto dimA(nc_read_blitz<int,1>(ncio.nc, "dim"+sgrids[1]));

            // ------------ wM
            auto wM_d(nc_read_blitz<double,1>(ncio.nc, BvA+".wM"));
            for (int i=0; i<wM_d.extent(0); ++i)
                wM.add({dimB(i)}, wM_d(i));


            // ------------ M
            std::array<long,2> shape;
            auto _dimB(ncio.nc->getVar("dim"+sgrids[0]));
            get_or_put_att(_dimB, 'r', "sparse_extent", "long", &shape[0], 1);
            auto _dimA(ncio.nc->getVar("dim"+sgrids[1]));
            get_or_put_att(_dimA, 'r', "sparse_extent", "long", &shape[1], 1);

            // Transfer full matrix shape meta-data
            M.set_shape(shape);
            wM.set_shape({shape[0]});
            Mw.set_shape({shape[1]});
            
//                std::array<long,2>{dimB.sparse_extent(), dimA.sparse_extent()});

            auto indices_d(nc_read_blitz<int,2>(ncio.nc, BvA+".M.indices"));
            auto values_d(nc_read_blitz<double,1>(ncio.nc, BvA+".M.values"));
            for (int i=0; i<values_d.extent(0); ++i) {
                int const iB_s = dimB(indices_d(i,0));
                int const iA_s = dimA(indices_d(i,1));
                if (scale) {
                    M.add({iB_s, iA_s}, values_d(i) * sM_s(iB_s));
                } else {
                    M.add({iB_s, iA_s}, values_d(i));
                }
            }

            /** Check that AvE is local */
            if (BvA == "AvE") {
                for (int i=0; i<values_d.extent(0); ++i) {
                    int const iAd = indices_d(i,0);
                    int const iA = dimB(iAd);
                    int const iEd = indices_d(i,1);
                    int const iE = dimA(iEd);

                    int const iA2 = iE % 12960;

                    if (iA != iA2) (*icebin_error)(-1,
                        "%d: AvE is not local!  iA=%d, iA2=%d, iE=%d\n",
                        iGlobal, iA, iA2, iE);
                    ++ iGlobal;

                }
            }

            // ------------------ Mw
            auto Mw_d(nc_read_blitz<double,1>(ncio.nc, BvA+".Mw"));
            for (int i=0; i<Mw_d.extent(0); ++i) Mw.add({dimA(i)}, Mw_d(i));
        }

        // Check
        if (ret.M.nnz() != nnz) (*icebin_error)(-1, "Bad count: %d vs %d", iM, nnz);
    }    // Finish off accumulators

    // Write it out
    printf("---- Writing output to %s\n", ofname.c_str());
    {NcIO ncio(ofname, ofmode);
        meta.ncio(ncio);
        ret.ncio(ncio, BvA);
    }
}

int main(int argc, char **argv)
{
    ParseArgs args(argc, argv);

    // Parse the output name
    const char *ifnames_cstr = args.ifnames[0].c_str();
    char const *dash = strrchr(ifnames_cstr, '-');
    std::string ofname(ifnames_cstr, dash);

    char ofmode = 'w';
    for (auto const &sgrids : args.matrix_specs) {
        combine_chunks(args.ifnames, ofname, ofmode, sgrids, args.scale);
        ofmode = 'a';
    }

    return 0;

}
