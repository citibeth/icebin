#include <cstdio>
#include <iostream>
#include <functional>
#include <icebin/GCMRegridder.hpp>
#include <icebin/IceRegridder_L0.hpp>
#include <icebin/Grid.hpp>
#include <spsparse/netcdf.hpp>

using namespace std;
using namespace netCDF;
using namespace ibmisc;
using namespace std::placeholders;  // for _1, _2, _3...
using namespace spsparse;

namespace icebin {

static double const nan = std::numeric_limits<double>::quiet_NaN();

// -----------------------------------------------------
IceRegridder::IceRegridder() : interp_style(InterpStyle::Z_INTERP), _name("icesheet") {}

IceRegridder::~IceRegridder() {}



// -------------------------------------------------------------
// ==============================================================
// Different vector spaces:
//      Description                    Domain
// ---------------------------------------------
// A  = Atmosphere grid                sphere
// Ap = Projected atmosphere grid      plane
// E  = Elevation grid                 sphere
// Ep = Projected elevation grid       plane
// I  = Ice grid                       plane
//
// Write out the parts that this class computed --- so we can test/check them
// -------------------------------------------------------------
/** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
NOTE: wAvAp == sApvA */
void IceRegridder::sApvA(MakeDenseEigenT::AccumT &&w) const
{
    for (int id=0; id < gcm->agridA.dim.dense_extent(); ++id) {
        auto index = gcm->agridA.dim.to_sparse(id);
        w.add({index, index}, gcm->agridA.native_area(id) / gridA_proj_area(id));
    }
}

/** Produces the diagonal matrix [Elevation projected] <-- [Elevation]
NOTE: wAvAp == sApvA */
void IceRegridder::sEpvE(MakeDenseEigenT::AccumT &&w) const
{
    for (int id=0; id < gcm->agridA.dim.dense_extent(); ++id) {
        auto index = gcm->agridA.dim.to_sparse(id);

        long nhc = gcm->nhc(index);
        long tuple[2] = {index, 0};
        long &ihp(tuple[1]);
        for (ihp=0; ihp<nhc; ++ihp) {
            long indexE = gcm->indexingHC.tuple_to_index(tuple);
            w.add({indexE, indexE}, gcm->agridA.native_area(id) / gridA_proj_area(id));
        }
    }
}
// -------------------------------------------------------------
#if 0
void IceRegridder::clear()
{
    agridI.clear();
    exgridI.clear();
}
#endif
// -------------------------------------------------------------
void IceRegridder::ncio(NcIO &ncio, std::string const &vname)
{
    if (ncio.rw == 'r') {
        agridI.ncio(ncio, vname + ".gridI");
        aexgrid.ncio(ncio, vname + ".exgrid");
    }

    auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
    get_or_put_att(info_v, ncio.rw, "name", _name);
    get_or_put_att_enum(info_v, ncio.rw, "interp_style", interp_style);

    ncio_blitz(ncio, gridA_proj_area, vname + ".gridA_proj_area", "double",
        get_or_add_dims(ncio, gridA_proj_area, {"gridA.ndata"}));
    agridI.ncio(ncio, vname + ".agridI");
    aexgrid.ncio(ncio, vname + ".aexgrid");

}

void IceRegridder::init(
    std::string const &name,
    AbbrGrid const &agridA,
    Grid const &fgridA,
    AbbrGrid const &&_agridI,
    AbbrGrid const &&_aexgrid,
    InterpStyle _interp_style)
{
    agridI = std::move(_agridI);    // convert Grid -> AbbrGrid
    aexgrid = std::move(_aexgrid);  // convert Grid -> AbbrGrid
    _name = (name != "" ? name : agridI.name);
    interp_style = _interp_style;

    if (agridI.sproj == "") {
        // No projection; projected and unproject area are the same
        gridA_proj_area.reference(agridA.native_area);
    } else {
        // Use a projection
        gridA_proj_area.reference(blitz::Array<double,1>(agridA.dim.dense_extent()));
        ibmisc::Proj_LL2XY proj(agridI.sproj);
        for (auto cell=fgridA.cells.begin(); cell != fgridA.cells.end(); ++cell) {
            int const is = cell->index;    // sparse index
            int const id = agridA.dim.to_dense(is);
            gridA_proj_area(id) = cell->proj_area(&proj);
        }
    }
}

std::unique_ptr<IceRegridder> new_ice_regridder(IceRegridder::Type type)
{
    switch(type.index()) {
        case IceRegridder::Type::L0 :
            return std::unique_ptr<IceRegridder>(new IceRegridder_L0);
        break;
        default :
            (*icebin_error)(-1,
                "Unknown IceRegridder::Type %s", type.str());
        break;
    }
}
std::unique_ptr<IceRegridder> new_ice_regridder(NcIO &ncio, std::string const &vname)
{
    std::string vn(vname + ".gridI.info");
    auto gridI_info_v = get_or_add_var(ncio, vn, "int64", {});

    GridParameterization parameterization;
    get_or_put_att_enum(gridI_info_v, ncio.rw, "parameterization", parameterization);
    return new_ice_regridder(parameterization);
}

// ---------------------------------------------------------------------
/** Made for binding... */
static bool in_good(std::unordered_set<long> const *set, long index_c)
{
    return (set->find(index_c) != set->end());
}

void IceRegridder::filter_cellsA(std::function<bool (long)> const &useA)
{

  // Figure out which cells to keep

    // List of cells in gridI / exgrid that overlap a cell we want to keep
    std::unordered_set<long> good_index_gridI;
    std::unordered_set<long> good_index_exgrid;


    std::unordered_set<int> good_j;
    for (int id=0; id<aexgrid.dim.dense_extent(); ++id) {
        auto const is = aexgrid.dim.to_sparse(id);
        if (useA(is)) {
            good_index_gridI.insert(aexgrid.ijk(id,1));    // j
            good_index_exgrid.insert(is);
        }
    }

    // Remove unneeded cells from gridI
    agridI.filter_cells(std::bind(&in_good, &good_index_gridI, _1));
    aexgrid.filter_cells(std::bind(&in_good, &good_index_exgrid, _1));
}
// ================================================================
// ==============================================================
extern void linterp_1d(
    std::vector<double> const &xpoints,
    double xx,
    int *indices, double *weights)  // Size-2 arrays
{
    int n = xpoints.size();


    // This is the point ABOVE our value.
    // (i0 = i1 - 1, xpoints[i0] < xx <= xpoints[i1])
    // See: http://www.cplusplus.com/reference/algorithm/lower_bound/
    int i1 = lower_bound(xpoints.begin(), xpoints.end(), xx) - xpoints.begin();

    if (i1 <= 0) i1 = 1;
    if (i1 >= n) i1 = n-1;

    int i0 = i1-1;
    indices[0] = i0;
    indices[1] = i1;
    double ratio = (xx - xpoints[i0]) / (xpoints[i1] - xpoints[i0]);

if (ratio < 0.0 || ratio > 1.0) {
    printf("BAD WEIGHTS: %g [%d]", xx, n);
    for (int i=0; i<n; ++i) printf(" %f", xpoints[i]);
    printf("\n");
}

    weights[0] = (1.0 - ratio);
    weights[1] = ratio;
}
// ========================================================================



} // namespace
