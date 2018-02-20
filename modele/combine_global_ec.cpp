#include <cstdio>
#include <ibmisc/netcdf.hpp>
#include <icebin/error.hpp>

using namespace std;
using namespace ibmisc;
using namespace netCDF;
using namespace icebin;

static double const NaN = std::numeric_limits<double>::quiet_NaN();

void combine_chunks(
    std::vector<std::string> const &ifnames,    // Names of input chunks
    std::string const &ofname,
    std::array<std::string,2> const &sgrids)    // {"B", "A"} --> matrix BvA
{
    printf("======== BEGIN combine_chunks(%sv%s)\n", sgrids[0].c_str(), sgrids[1].c_str());

    // Names of the variables we will read/write
    std::string BvA = sgrids[0] + "v" + sgrids[1];

    // Get total sizes
    int nnz = 0;        // = number non-zero]
    std::array<size_t,2> sparse_extents;
    for (std::string const &ifname : ifnames) {
        NcIO ncio(ifname, 'r');
        nnz += ncio.nc->getDim(BvA+".M.size").getSize();

        for (int i=0; i<2; ++i) {
            NcVar nc_dimA(ncio.nc->getVar("dim" + sgrids[i]));
            nc_dimA.getAtt("sparse_extent").getValues(&sparse_extents[i]);
        }
    }

    // Allocate
    blitz::Array<int,2> indices(nnz,2);
    indices = -1;
    blitz::Array<double,1> wM(sparse_extents[0]);
    wM = NaN;
    blitz::Array<double,1> values(nnz,1);
    values = NaN;
    blitz::Array<double,1> Mw(sparse_extents[1]);
    Mw = NaN;

    // Agreate them together
    int iM=0;
    for (std::string const &ifname : ifnames) {
        NcIO ncio(ifname, 'r');

        auto dimB(nc_read_blitz<int,1>(ncio.nc, "dim"+sgrids[0]));
        auto dimA(nc_read_blitz<int,1>(ncio.nc, "dim"+sgrids[1]));

        auto wM_d(nc_read_blitz<double,1>(ncio.nc, BvA+".wM"));
        for (int i=0; i<wM_d.extent(0); ++i) wM(dimB(i)) = wM_d(i);

        auto indices_d(nc_read_blitz<int,1>(ncio.nc, BvA+".M.indices"));
        auto values_d(nc_read_blitz<double,1>(ncio.nc, BvA+".M.values"));
        for (int i=0; i<values.extent(0); ++i) {
            values(iM) = values_d(i);
            indices(iM,0) = indices_d(i,0);
            indices(iM,1) = indices_d(i,1);
            ++iM;
        }

        auto Mw_d(nc_read_blitz<double,1>(ncio.nc, BvA+".Mw"));
        for (int i=0; i<Mw_d.extent(0); ++i) Mw(dimA(i)) = Mw_d(i);
    }

    // Check
    if (iM != nnz) (*icebin_error)(-1, "Bad count: %d vs %d", iM, nnz);

    // Write it out
    {NcIO ncio(ofname, 'a');
        // Dimensions
        auto dims(get_or_add_dims(ncio,
            {"dim"+sgrids[0], "dim"+sgrids[1]},
            {sparse_extents[0], sparse_extents[1]}));
        auto nnz_d(get_or_add_dim(ncio, BvA+".nnz", nnz));

        // Variables
        ncio_blitz(ncio, wM, BvA+".wM", "double", {dims[0]});
        ncio_blitz(ncio, indices, BvA+".indices", "int", {nnz_d});
        ncio_blitz(ncio, values, BvA+".values", "double", {nnz_d});
        ncio_blitz(ncio, Mw, BvA+".Mw", "double", {dims[1]});

        ncio.flush();
    }

}

int main(int argc, char **argv)
{
    vector<string> ifnames {"global_ec.nc-00", "global_ec.nc-01", "global_ec.nc-02", "global_ec.nc-03", "global_ec.nc-04"};
    string ofname("global_ec.nc");
    vector<array<string,2>> sgridss {
        {"I","E"}, {"E","I"},
        {"I","A"}, {"A","I"},
        {"A","E"}, {"E","A"}};

    for (auto const &sgrids : sgridss) {
        combine_chunks(ifnames, ofname, sgrids);
    }
}
