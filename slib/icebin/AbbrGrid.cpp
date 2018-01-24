#include <icebin/AbbrGrid.hpp>
#include <icebin/Grid.hpp>
#include <ibmisc/netcdf.hpp>

using namespace ibmisc;

namespace icebin {

void AbbrGrid::operator=(Grid const &g)
{
    spec = g.spec->clone();
    coordinates = g.coordinates;
    parameterization = g.parameterization;
    indexing = g.indexing;
    name = g.name;
    sproj = g.sproj;

    // Allocate
    auto nd = g.ndata();
    ijk.reference(blitz::Array<int,2>(nd,3));
    native_area.reference(blitz::Array<double,1>(nd));
    centroid.reference(blitz::Array<double,2>(nd,2));

    // Copy info into AbbrGrid
    int id=0;
    ibmisc::Proj_LL2XY proj(g.sproj);
    for (auto cell=g.cells.begin(); cell != g.cells.end(); ++cell) {
        dim.add_dense(cell->index);
        ijk[0](id) = cell->i;
        ijk[1](id) = cell->j;
        ijk[2](id) = cell->k;
        native_area(id) = cell->native_area;
        auto ctr(cell->centroid());
        centroid(id,0) = ctr.x;
        centroid(id,1) = ctr.y;
    }
}

void AbbrGrid::filter_cells(std::function<bool(long)> const &keep_fn)
{

    SparseSet<long,int> &dim0(dim);
    SparseSet<long,int> dim1(dim0.sparse_extent());

    // Determine new dimension set
    blitz::Array<double,1> dim0v1(dim0.dense_extent());    // Convert dense dim0 <-- dense dim1 indices
    for (int id0=0; id0<dim0.dense_extent(); ++id0) {
        int is = dim0.to_sparse(id);
        if (keep_fn(is)) {
            auto id1 = dim1.add_dense(is);
            dim0v1(id1) = id0;
        }
    }

    // Keep old arrays for now
    blitz::Array<int,2> ijk0(ijk);
    blitz::Array<double,1> native_area0(native_area);
    blitz::Array<double,2> centroid0(centroid);

    // Allocate new arrays
    auto const N = dim1.dense_extent();
    ijk.reference(blitz::Array<int,2>(N,3));
    native_area.reference(blitz::Array<double,1>(N));
    centroid.reference(blitz::Array<double,2>(N,2));

    // Copy over
    for (int id1=0; id1<dim1.dense_extent(); ++id1) {
        auto const id0 = dim0v1(id1);
        ijk(id1,0) = ijk0(id0,0);
        ijk(id1,1) = ijk0(id0,1);
        ijk(id1,2) = ijk0(id0,2);
        native_area(id1) = native_area0(id0);
        centroid(id1,0) = centroid0(id0,0);
        centroid(id1,1) = centroid0(id0,1);
    }

    // Put new dim into place
    dim = std::move(dim1);

    return agrid1;
}


void AbbrGrid::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    ncio_grid_spec(ncio, spec, vname);

    auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});
    get_or_put_att_enum(info_v, ncio.rw, "type", type);
    get_or_put_att_enum(info_v, ncio.rw, "coordinates", coordinates);
    get_or_put_att_enum(info_v, ncio.rw, "parameterization", parameterization);
    get_or_put_att(info_v, ncio.rw, "name", name);
    get_or_put_att(info_v, ncio.rw, "sproj", sproj);

    indexing.ncio(ncio, vname + ".indexing");
    dim.ncio(ncio, vname + ".dim");

    auto sparse_extent_d(get_or_add_dim(ncio,
        vname+".dim.sparse_extent", ijk.extent(0)));
    auto three_d(get_or_add_dim(ncio, "three", 3));
    auto two_d(get_or_add_dim(ncio, "two", 2));

    ncio_blitz_alloc(ncio, ijk, vname + ".ijk", "int",
        {sparse_extent_d, three_d});
    ncio_blitz_alloc(ncio, native_area, vname + ".native_area", "double",
        {sparse_extent_d});
    ncio_blitz_alloc(ncio, centroid, vname + ".ijk", "double",
        {sparse_extent_d, two_d});

}


}    // namespace
