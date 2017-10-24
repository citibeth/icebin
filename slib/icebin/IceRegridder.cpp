#include <cstdio>
#include <iostream>
#include <functional>
#include <icebin/GCMRegridder.hpp>
#include <icebin/IceRegridder_L0.hpp>
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
void IceRegridder::sApvA(MakeDenseEigenT::AccumT &w) const
{
    ibmisc::Proj_LL2XY proj(gridI->sproj);
    for (auto cell=gcm->gridA->cells.begin(); cell != gcm->gridA->cells.end(); ++cell) {
        w.add({cell->index, cell->index}, cell->native_area / cell->proj_area(&proj));
    }
}

/** Produces the diagonal matrix [Elevation projected] <-- [Elevation]
NOTE: wAvAp == sApvA */
void IceRegridder::sEpvE(MakeDenseEigenT::AccumT &w) const
{
    ibmisc::Proj_LL2XY proj(gridI->sproj);
    for (auto cell=gcm->gridA->cells.begin(); cell != gcm->gridA->cells.end(); ++cell) {
        long nhc = gcm->nhc(cell->index);
        long tuple[2] = {cell->index, 0};
        long &ihp(tuple[1]);
        for (ihp=0; ihp<nhc; ++ihp) {
            long indexE = gcm->indexingHC.tuple_to_index(tuple);
            w.add({indexE, indexE}, cell->native_area / cell->proj_area(&proj));
        }
    }
}
// -------------------------------------------------------------
void IceRegridder::clear()
{
    gridI.reset();
    exgrid.reset();
}
// -------------------------------------------------------------
void IceRegridder::ncio(NcIO &ncio, std::string const &vname, bool rw_full)
{
printf("BEGIN IceRegridder::ncio(%s, %d)\n", vname.c_str(), rw_full);
    if (ncio.rw == 'r') {
        clear();
        gridI = new_grid(ncio, vname + ".gridI");
        exgrid = new_grid(ncio, vname + ".exgrid");
    }

    auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});
    get_or_put_att(info_v, ncio.rw, "name", _name);
    get_or_put_att_enum(info_v, ncio.rw, "interp_style", interp_style);

    gridI->ncio(ncio, vname + ".gridI", rw_full);
    exgrid->ncio(ncio, vname + ".exgrid", rw_full);

printf("END IceRegridder::ncio(%s, %d)\n", vname.c_str(), rw_full);
}

void IceRegridder::init(
    std::string const &name,
    std::unique_ptr<Grid> &&_gridI,
    std::unique_ptr<Grid> &&_exgrid,
    InterpStyle _interp_style)
{
    _name = (name != "" ? name : gridI->name);
    gridI = std::move(_gridI);
    exgrid = std::move(_exgrid);
    interp_style = _interp_style;
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

    Grid::Parameterization parameterization;
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
    for (auto excell = exgrid->cells.begin(); excell != exgrid->cells.end(); ++excell) {
        int index1 = excell->i;
        if (useA(index1)) {
            good_index_gridI.insert(excell->j);
            good_index_exgrid.insert(excell->index);
        }
    }

    // Remove unneeded cells from gridI
    gridI->filter_cells(std::bind(&in_good, &good_index_gridI, _1));
    exgrid->filter_cells(std::bind(&in_good, &good_index_exgrid, _1));
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
