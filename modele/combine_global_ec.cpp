#include <cstdio>
#include <ibmisc/netcdf.hpp>
#include <icebin/error.hpp>
#include <prettyprint.hpp>

using namespace std;
using namespace ibmisc;
using namespace netCDF;
using namespace icebin;

static double const NaN = std::numeric_limits<double>::quiet_NaN();


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



template<int RANK>
void _write_reshaped(
    NcIO &ncio,
    blitz::Array<double,1> const &var,
    std::string const &vname,
    std::string const &stype,
    std::vector<int> const &shape,
    std::string const &sdim,
    std::array<std::string,RANK> const &sdim_exts)
{
    blitz::TinyVector<int,RANK> bshape();

    auto _var2(reshape<double,1,RANK>(var, vector_to_tiny<int,RANK>(shape)));
    blitz::Array<double,RANK> &var2(ncio.tmp.copy(_var2));    // Copy constructor, just reference

printf("_write_reshaped RANK=%d\n", RANK);
printf("var = [%g %g %g...]\n", var(0), var(1), var(2));


    std::vector<std::string> dim_names;
    for (int i=0; i<RANK; ++i)
        dim_names.push_back(sdim + sdim_exts[i]);
    auto dims(get_or_add_dims(ncio, var2, dim_names));
    ncio_blitz(ncio, var2, vname, stype, dims);
}


void write_reshaped(
    NcIO &ncio,
    blitz::Array<double,1> const &var,
    std::string const &vname,
    std::string const &stype,
    std::vector<int> const &shape,
    int sparse_extent,
    std::string const &sdim)
{
        switch(shape.size()) {
            case 2:
                _write_reshaped<2>(ncio, var, vname, "double",
                    shape, sdim, {".jm", ".im"});
            break;
            case 3:
                _write_reshaped<3>(ncio, var, vname, "double",
                    shape, sdim, {".nhc", ".jm", ".im"});
            break;
            default:
                _write_reshaped<1>(ncio, var, vname, "double",
                    {sparse_extent}, sdim, {""});
        }
}



void combine_chunks(
    std::vector<std::string> const &ifnames,    // Names of input chunks
    std::string const &ofname,
    char ofmode,    // 'w' or 'a' for write mode of output file
    std::array<std::string,2> const &sgrids)    // {"B", "A"} --> matrix BvA
{
    printf("======== BEGIN combine_chunks(%sv%s)\n", sgrids[0].c_str(), sgrids[1].c_str());

    // Names of the variables we will read/write
    std::string BvA = sgrids[0] + "v" + sgrids[1];

    // Get total sizes
    int nnz = 0;        // = number non-zero]
    std::array<size_t,2> sparse_extents;
    std::array<std::vector<int>,2> shapes;
    std::vector<size_t> sizes;    // For printing
    for (std::string const &ifname : ifnames) {
        NcIO ncio(ifname, 'r');
        size_t sz = ncio.nc->getDim(BvA+".M.nnz").getSize();
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


    // Allocate
    blitz::Array<int,2> indices(nnz,2);
    indices = -1;
    blitz::Array<double,1> wM(sparse_extents[0]);
    wM = NaN;
    blitz::Array<double,1> values(nnz);
    values = NaN;
    blitz::Array<double,1> Mw(sparse_extents[1]);
    Mw = NaN;

    // Agreate them together
    int iM=0;
    for (std::string const &ifname : ifnames) {
        printf("--- Reading %s\n", ifname.c_str());
        NcIO ncio(ifname, 'r');

        auto dimB(nc_read_blitz<int,1>(ncio.nc, "dim"+sgrids[0]));
        auto dimA(nc_read_blitz<int,1>(ncio.nc, "dim"+sgrids[1]));

        auto wM_d(nc_read_blitz<double,1>(ncio.nc, BvA+".wM"));
        for (int i=0; i<wM_d.extent(0); ++i) wM(dimB(i)) = wM_d(i);

        auto indices_d(nc_read_blitz<int,2>(ncio.nc, BvA+".M.indices"));
        auto values_d(nc_read_blitz<double,1>(ncio.nc, BvA+".M.values"));
        for (int i=0; i<values_d.extent(0); ++i) {
            if (iM >= nnz) (*icebin_error)(-1,
                "iM=%d is too large (nnz=%d)\n", iM, nnz);
            values(iM) = values_d(i);
            indices(iM,0) = dimB(indices_d(i,0));
            indices(iM,1) = dimA(indices_d(i,1));
            ++iM;
        }

        auto Mw_d(nc_read_blitz<double,1>(ncio.nc, BvA+".Mw"));
        for (int i=0; i<Mw_d.extent(0); ++i) Mw(dimA(i)) = Mw_d(i);
    }

    // Check
    if (iM != nnz) (*icebin_error)(-1, "Bad count: %d vs %d", iM, nnz);

    // Write it out
    printf("---- Writing output to %s\n", ofname.c_str());
    {NcIO ncio(ofname, ofmode);
        // Dimensions
        auto nnz_d(get_or_add_dim(ncio, BvA+".nnz", nnz));
        auto two_d(get_or_add_dim(ncio, "two", 2));


        // Variables
        write_reshaped(ncio, wM, "w"+BvA, "double", shapes[0], sparse_extents[0], "dim"+sgrids[0]);
        ncio_blitz(ncio, indices, BvA+".indices", "int", {nnz_d, two_d});
        ncio_blitz(ncio, values, BvA+".values", "double", {nnz_d});
        write_reshaped(ncio, Mw, BvA+"w", "double", shapes[1], sparse_extents[1], "dim"+sgrids[1]);

        ncio.flush();
    }

}

int main(int argc, char **argv)
{
    // Save args as C++ vector
    vector<string> ifnames;
    for (int i=1; i<argc; ++i) ifnames.push_back(string(argv[i]));

    // Parse the output name
    char *dash = rindex(argv[1], '-');
    std::string ofname(argv[1], dash);

    vector<array<string,2>> sgridss {
        {"I","E"}, {"E","I"},
        {"I","A"}, {"A","I"},
        {"A","E"}};

    char ofmode = 'w';
    for (auto const &sgrids : sgridss) {
        combine_chunks(ifnames, ofname, ofmode, sgrids);
        ofmode = 'a';
    }

    return 0;
}
