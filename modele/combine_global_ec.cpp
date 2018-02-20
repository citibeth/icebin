#include <cstdio>
#include <ibmisc/netcdf.hpp>

using namespace ibmisc;

void combine_chunks(
    std::vector<std::string> const &ifnames,    // Names of input chunks
    std::vector<std::string> const &ofname
    std::array<std::string> const &sgrids)    {"B", "A"} --> matrix BvA
{
    printf("======== BEGIN combine_chunks(%sv%s)\n", sgrids[0], sgrids[1]);

    // Names of the variables we will read/write
    std::string vn_dimB("dim" + Bgrid);
    std::string vn_dimA("dim" + Agrid);
    std::string vn_BvA = Bgrid + "v" + Agrid;

#if 0
    std::string vn_wM(strprintf("%sv%s.wM", Bgrid.c_str(), Agrid.c_str()));
    std::string vn_M(strprintf("%sv%s.M", Bgrid.c_str(), Agrid.c_str()));
    std::string vn_Mw(strprintf("%sv%s.Mw", Bgrid.c_str(), Agrid.c_str()));
    std::string vn_info(strprintf("%sv%s.M.info", Bgrid.c_str(), Agrid.c_str()));
#endif

    // Get total sizes
    int nnz = 0;        // = number non-zero]
    std::array<size_t,2> sparse_extents;
    for (std::string const &ifname : ifnames) {
        NcIO ncio(ifname, 'r');
        nnz += ncio.nc->getDim(vn_BvA + ".M.size").getSize();

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

        auto wM_d(nc_read_blitz<double,1>(ncio.nc, vn+bva + ".wM"));
        for (int i=0; i<wM_d.extent(0); ++i) wM(dimB(i)) = wM_d(i);

        auto indices_d(nc_read_blitz<int,1>(ncio.ncc, vn+bva + ".M.indices"));
        auto values_d(nc_read_blitz<double,1>(ncio.nc, vn+bva + ".M.values"));
        for (int i=0; i<values.extent(0); ++i) {
            values(iM) = values_d(i);
            indices(iM,0) = indices_d(i,0);
            indices(iM,1) = indices_d(i,1);
            ++iM;
        }

        auto Mw_d(nc_read_blitz<double,1>(ncio.nc, vn+bva + ".Mw"));
        for (int i=0; i<Mw_d.extent(0); ++i) Mw(dimA(i)) = Mw_d(i);
    }

    // Check
    if (iM != nnz) (*icebin_error)(-1, "Bad count: %d vs %d", iM, nnz);

    // Write it out
    {NcIO ncio(ofname, 'a');
        // Dimensions
        auto dims(get_or_add_dims(ncio,
            {"dim"+sgrids[0], "dim"+sgrids[1]},
            sparse_extents));
        auto nnz_d(get_or_add_dim(ncio, vn_BvA+".nnz", nnz));

        // Variables
        ncio_blitz(ncio, wM, vn_BvA+".wM", "double", {dims[0]});
        ncio_blitz(ncio, indices, vn_BvA+".indices", "int", {nnz_d});
        ncio_blitz(ncio, values, vn_BvA+".values", "double", {nnz_d});
        ncio_blitz(ncio, Mw, vn_BvA+".Mw", "double", {dims[1]});

        ncio.flush();
    }

}

int main(int argc, char **argv)
{
    std::vector<std::string> ifnames {"global_ec.nc-00", "global_ec.nc-01", "global_ec.nc-02", "global_ec.nc-03", "global_ec.nc-04"};
    std::string ofname("global_ec.nc");
    std::vector<std::array<std::string>> sgridss {
        {"I","E"}, {"E","I"},
        {"I","A"}, {"A","I"},
        {"A","E"}, {"E","A"}};

    for (auto const &sgrids : sgridss) {
        combine_chunks(ifnames, ofnames, sgrids);
    }
}
