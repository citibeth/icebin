#include <icebin/AbbrGrid.hpp>
#include <icebin/Grid.hpp>
#include <ibmisc/netcdf.hpp>

using namespace ibmisc;
using namespace spsparse;

namespace icebin {

AbbrGrid::AbbrGrid(
    std::unique_ptr<GridSpec> &&_spec,
    GridCoordinates _coordinates,
    GridParameterization _parameterization,
    ibmisc::Indexing _indexing,
    std::string const &_name,
    std::string const &_sproj,
    spsparse::SparseSet<long,int> &&_dim,
    blitz::Array<int,2> const &_ijk,
    blitz::Array<double,1> const &_native_area,
    blitz::Array<double,2> const &_centroid_xy)
: spec(std::move(_spec)),
coordinates(_coordinates), parameterization(_parameterization),
indexing(_indexing), name(_name), sproj(_sproj),
dim(std::move(_dim)),
ijk(_ijk), native_area(_native_area), centroid_xy(_centroid_xy)
{}


void AbbrGrid::clear()
{
    spec.reset();
    name = "";
    sproj = "";
    dim.clear();
    ijk.free();
    native_area.free();
    centroid_xy.free();
}


/** Convert from Grid */
AbbrGrid::AbbrGrid(Grid const &g) :
    spec(g.spec),
    coordinates(g.coordinates),
    parameterization(g.parameterization),
    indexing(g.indexing),
    name(g.name),
    sproj(g.sproj),
    dim(g.ndata())    // sparse_extent
{
    // Allocate
    auto nd = g.nrealized();    // dense extent
    ijk.reference(blitz::Array<int,2>(nd,3));
    native_area.reference(blitz::Array<double,1>(nd));
    if (g.coordinates == GridCoordinates::XY) {
        centroid_xy.reference(blitz::Array<double,2>(nd,2));
    }

    // Copy info into AbbrGrid
    ibmisc::Proj_LL2XY proj(g.sproj);
    for (auto cell=g.cells.begin(); cell != g.cells.end(); ++cell) {
        int id = dim.add_dense(cell->index);    // Dense index
        if (id >= nd) (*icebin_error)(-1,
            "Index out of range: %d vs %d", id, nd);

        ijk[0](id) = cell->i;
        ijk[1](id) = cell->j;
        ijk[2](id) = cell->k;
        native_area(id) = cell->native_area;
        if (g.coordinates == GridCoordinates::XY) {
            auto ctr(cell->centroid());
            centroid_xy(id,0) = ctr.x;
            centroid_xy(id,1) = ctr.y;
        }
    }
}

void AbbrGrid::filter_cells(std::function<bool(long)> const &keep_fn)
{

    SparseSet<long,int> &dim0(dim);
    SparseSet<long,int> dim1(dim0.sparse_extent());

    // Determine new dimension set
    blitz::Array<int,1> dim0v1(dim0.dense_extent());    // Convert dense dim0 <-- dense dim1 indices
    for (int id0=0; id0<dim0.dense_extent(); ++id0) {
        int is = dim0.to_sparse(id0);
        if (keep_fn(is)) {
            auto id1 = dim1.add_dense(is);
            dim0v1(id1) = id0;
        }
    }

    // Keep old arrays for now
    blitz::Array<int,2> ijk0(ijk);
    blitz::Array<double,1> native_area0(native_area);
    blitz::Array<double,2> centroid_xy0(centroid_xy);

    // Allocate new arrays
    auto const N = dim1.dense_extent();
    ijk.reference(blitz::Array<int,2>(N,3));
    native_area.reference(blitz::Array<double,1>(N));
    centroid_xy.reference(blitz::Array<double,2>(N,2));

    // Copy over
    for (int id1=0; id1<dim1.dense_extent(); ++id1) {
        auto const id0 = dim0v1(id1);
        ijk(id1,0) = ijk0(id0,0);
        ijk(id1,1) = ijk0(id0,1);
        ijk(id1,2) = ijk0(id0,2);
        native_area(id1) = native_area0(id0);
        centroid_xy(id1,0) = centroid_xy0(id0,0);
        centroid_xy(id1,1) = centroid_xy0(id0,1);
    }

    // Put new dim into place
    dim = std::move(dim1);
}


void AbbrGrid::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    ncio_grid_spec(ncio, spec, vname);

    auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
    get_or_put_att_enum(info_v, ncio.rw, "coordinates", coordinates);
    get_or_put_att_enum(info_v, ncio.rw, "parameterization", parameterization);
    get_or_put_att(info_v, ncio.rw, "name", name);
    get_or_put_att(info_v, ncio.rw, "sproj", sproj);

    indexing.ncio(ncio, vname + ".indexing");

    // Store dim; retrieve dimension from it
    dim.ncio(ncio, vname + ".dim");

    auto dense_extent_d(get_or_add_dim(ncio,
        vname+".dim.dense_extent", ijk.extent(0)));    // extent ignored on read
    auto three_d(get_or_add_dim(ncio, "three", 3));
    auto two_d(get_or_add_dim(ncio, "two", 2));

    ncio_blitz_alloc(ncio, ijk, vname + ".ijk", "int",
        {dense_extent_d, three_d});
    ncio_blitz_alloc(ncio, native_area, vname + ".native_area", "double",
        {dense_extent_d});
    ncio_blitz_alloc(ncio, centroid_xy, vname + ".centroid_xy", "double",
        {dense_extent_d, two_d});

}

// ==============================================================================
// We only need to define these because blitz::Array does not follow STL conventions
// and is not movable.
void AbbrGrid::operator=(AbbrGrid &&other)
{
    // Standard stuff
    spec = std::move(other.spec);
    coordinates = other.coordinates;
    parameterization = other.parameterization;
    indexing = std::move(other.indexing);
    name = std::move(other.name);
    sproj = std::move(other.sproj);
    dim = std::move(other.dim);

    // Non-standard stuff, why we have to write our own operator=(&&)
    // We don't really move here, just copy the shared pointer.
    ijk.reference(other.ijk);
    native_area.reference(other.native_area);
    centroid_xy.reference(other.centroid_xy);
}

AbbrGrid::AbbrGrid(AbbrGrid &&other)
    { *this = std::move(other); }

void AbbrGrid::operator=(AbbrGrid const &other)
{
    // Standard stuff
    spec = other.spec;
    coordinates = other.coordinates;
    parameterization = other.parameterization;
    indexing = other.indexing;
    name = other.name;
    sproj = other.sproj;
    dim = other.dim;

    // Non-standard stuff, why we have to write our own operator=(&&)
    // We don't really move here, just copy the shared pointer.
    ijk.reference(other.ijk);
    native_area.reference(other.native_area);
    centroid_xy.reference(other.centroid_xy);
}

AbbrGrid::AbbrGrid(AbbrGrid const &other)
    { *this = other; }




}    // namespace
