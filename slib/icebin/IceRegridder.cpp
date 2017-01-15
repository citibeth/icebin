#include <cstdio>
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
void IceRegridder::sApvA(SparseTriplets<SparseMatrix> &w)
{
    ibmisc::Proj_LL2XY proj(gridI->sproj);
    for (auto cell=gcm->gridA->cells.begin(); cell != gcm->gridA->cells.end(); ++cell) {
        w.add({cell->index, cell->index}, cell->native_area / cell->proj_area(&proj));
    }
}

/** Produces the diagonal matrix [Elevation projected] <-- [Elevation]
NOTE: wAvAp == sApvA */
void IceRegridder::sEpvE(SparseTriplets<SparseMatrix> &w)
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
    elevI = nan;
}
// -------------------------------------------------------------
void IceRegridder::ncio(NcIO &ncio, std::string const &vname)
{
    if (ncio.rw == 'r') {
        clear();
        gridI = new_grid(ncio, vname + ".gridI");
        exgrid = new_grid(ncio, vname + ".exgrid");
    }

    auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});
    get_or_put_att(info_v, ncio.rw, "name", _name);
    get_or_put_att_enum(info_v, ncio.rw, "interp_style", interp_style);

    gridI->ncio(ncio, vname + ".gridI");
    exgrid->ncio(ncio, vname + ".exgrid");
    ncio_blitz(ncio, elevI, true, vname + ".elevI", "double",
        get_dims(ncio ,{vname + ".gridI.cells.nfull"}));
}

void IceRegridder::init(
    std::string const &name,
    std::unique_ptr<Grid> &&_gridI,
    std::unique_ptr<Grid> &&_exgrid,
    InterpStyle _interp_style,
    blitz::Array<double,1> &_elevI)
{
    _name = (name != "" ? name : gridI->name);
    gridI = std::move(_gridI);
    exgrid = std::move(_exgrid);
    interp_style = _interp_style;
    elevI.reference(_elevI);
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

/** Used to generate Ur matrices for Atm or Elevation grids.
This allows us to exploit an algebraic symmetry between A and E. */
struct UrAE {
    const long nfull;

    const std::function<void(MakeDenseEigenT::AccumT &)> GvAp;
    const std::function<void(MakeDenseEigenT::AccumT &)> sApvA;

    UrAE(long _nfull, GvAp_fn _GvAp, GvAp_fn _sApvA) : nfull(_nfull), GvAp(_GvAp), sApvA(_sApvA) {}
};



// ------------------------------------------------------------
std::function<bool(long)> in_sparse_fn(SparseSet const &dim)
    { return std::bind(&SparseSet::in_sparse, &dim, _1); }
// ------------------------------------------------------------
static std::unique_ptr<WeightedSparse> compute_AEvI(IceRegridder *regridder, bool scale, bool correctA, UrAE const &AE)
{
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse);
    SparseSet &dimA(ret->dims[0]);
    SparseSet &dimI(ret->dims[1]);
    SparseSet dimG;

    dimA.set_sparse_extent(AE.nfull);
    dimI.set_sparse_extent(regridder->nI());
    dimG.set_sparse_extent(regridder->nG());

    // ----- Get the Ur matrices (which determines our dense dimensions)
    MakeDenseEigenT GvI_m(
        std::bind(&IceRegridder::GvI, regridder, _1),
        SparsifyTransform::ADD_DENSE
        {&dimG, &dimI}, '.');
    MakeDenseEigenT ApvG_m(
        AE.GvAP,
        SparsifyTransform::ADD_DENSE,
        {&dimG, &dimA}, 'T');

    // ----- Convert to Eigen and multiply
    auto ApvG(ApvG_m.to_eigen());
    auto GvI(GvI_m.to_eigen());
    auto sGvI(sum_to_diagonal(GvI, 0, '-'));

    std::unique_ptr<EigenSparseMatrix> ApvI(
        new EigenSparseMatrix(ApvG * sGvI * GvI));

    // ----- Apply final scaling, and convert back to sparse dimension
    if (correctA) {
        // ----- Compute the final weight matrix
        auto wAvAp(MakeDenseEigenT(                   // diagonal
            AE.sApvA,
            SparseTransform::TO_DENSE_IGNORE_MISSING,
            {&dimA, &dimA}, '.').to_eigen());
        auto wApvI(sum_to_diagonal(*ApvI, 0, '+'));        // diagonal
        EigenSparseMatrix wAvI(wAvAp_e * wApvI_e);    // diagonal...

        // Sums rows of an Eigen matrix into a dense blitz::Array
        ret->weight.reference(sum(wAvI, 0, '+'));

        // Compute the main matrix
        if (scale) {
            // Get two diagonal Eigen scale matrices
            auto sAvAp(sum_to_diagonal(wAvAp, 0, '-'));
            auto sApvI(sum_to_diagonal(wApvI, 0, '-'));

            ret->M.reset(new EigenSparseMatrix(
                sAvAp * sApvI * *ApvI));    // AvI_scaled
        } else {
            ret->M = std::move(ApvI_e);
        }

    } else {

        // ----- Compute the final weight matrix
        auto wApvI_b(sum(*ApvI, 0, '+'));
        ret->weight.reference(wApvI_b);

        if (scale) {
            // Get two diagonal Eigen scale matrices
            auto sApvI(diag_matrix(wApvI_b, '-'));

            ret->M.reset(new EigenSparseMatrix(
                sApvI * *ApvI));    // ApvI_scaled
        } else {
            ret->M = std::move(ApvI);
        }
    }

    return ret;
}
// ---------------------------------------------------------

static std::unique_ptr<WeightedSparse> compute_IvAE(IceRegridder *regridder, bool scale, bool correctA, UrAE const &AE)
{
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse);
    auto &dimA(ret->dims[1]);
    auto &dimI(ret->dims[0]);
    SparseSet dimG;

    dimA.set_sparse_extent(AE.nfull);
    dimI.set_sparse_extent(regridder->nI());
    dimG.set_sparse_extent(regridder->nG());

    // ----- Get the Ur matrices (which determines our dense dimensions)
    MakeDenseEigenT GvAp_m(
        AE.GvAp,
        SparsifyTransform::ADD_DENSE,
        {&dimG, &dimA}, '.');
    MakeDenseEigenT IvG_m(
        std::bind(&IceRegridder::GvI, regridder, _1),
        SparsifyTransform::ADD_DENSE
        {&dimG, &dimI}, 'T');

    // ----- Convert to Eigen
    auto GvAp(GvAp_m.to_eigen());
    auto IvG(IvG_m.to_eigen());
    auto sGvAp(sum_to_diagonal(GvAp, 0, '-'));

    // Unscaled matrix
    std::unique_ptr<EigenSparseMatrixT> IvAp_e(
        new EigenSparseMatrixT(IvG * sGvAp * GvAp));

    // Get weight vector from IvAp_e
//    ret->weight.reference(sum(IvG, 0));  // testing only
    ret->weight.reference(sum(*IvAp, 0, '+'));


    // ----- Apply final scaling, and convert back to sparse dimension
    if (correctA) {
        // Scaling matrix
        auto sApvA(MakeDenseEigenT(
            AE.sApvA,
            SparseTransform::TO_DENSE_IGNORE_MISSING,
            {&dimA, &dimA}, '.').to_eigen());

        if (scale) {
            auto sIvAp(sum_to_diagonal(IvAp, 0, '-'));
            ret->M.reset(new EigenSparseMatrixT(sIvAp * *IvAp * sApvA));
        } else {
            ret->M.reset(new EigenSparseMatrixT(*IvAp * sApvA));
        }
    } else {
        if (scale) {
            auto sIvAp(scale_matrix(*IvAp, 0));
            ret->M.reset(new EigenSparseMatrixT(sIvAp * *IvAp));
        } else {
            ret->M = std::move(IvAp);
        }
    }

    return ret;
}

static std::unique_ptr<WeightedSparse> compute_EvA(IceRegridder *regridder, bool scale, bool correctA, UrAE const &E, UrAE const &A)
{
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse);
    auto &dimE(ret->dims[0]);
    auto &dimA(ret->dims[1]);
    SparseSet dimG;

    // ----- Get the Ur matrices (which determines our dense dimensions)

    MakeDenseEigenT GvAp_m(
        A.GvAp,
        SparsifyTransform::ADD_DENSE,
        {&dimG, &dimA}, '.');
    MakeDenseEigenT EpvG_m(
        E.GvAp,
        SparsifyTransform::ADD_DENSE,
        {&dimG, &dimE}, 'T');

    // ----- Convert to Eigen and multiply
    auto GvAp(GvAp_m.to_eigen());
    auto EpvG(EpvG_m.to_eigen());
    auto sGvAp(sum_to_diagonal(GvAp, 0, '-'));

    // Unweighted matrix
    std::unique_ptr<EigenSparseMatrix> EpvAp(
        new EigenSparseMatrix(EpvG * sGvAp * GvAp));
    auto wEpvAp_b(sum(*EpvAp,0,'+'));

    // ----- Apply final scaling, and convert back to sparse dimension
    if (correctA) {
        auto sApvA(MakeDenseEigenT(
            A.sApvA,
            SparseTransform::TO_DENSE_IGNORE_MISSING,
            {&dimA, &dimA}, '.').to_eigen());

        if (scale) {
            auto sEpvE(MakeDenseEigenT(
                E.sApvA,
                SparseTransform::TO_DENSE_IGNORE_MISSING,
                {&dimE, &dimE}, '.').to_eigen());


            auto sEvEp(sum_to_diagonal(sEpvE,0,'-'));
            auto sEpvAp(diag_matrix(wEpvAp_b,'-'));

            ret->M.reset(new EigenSparseMatrix(
                sEvEp * sEpvAp * *EpvAp * sApvA));    // EvA
        } else {
            ret->M.reset(new EigenSparseMatrix(*EpvAp * sApvA));
        }

        // ----- Compute the final weight (diagonal) matrix
        auto wEpvAp_b(sum(*EpvAp,0,'+'));
        ret->weight.reference(wEpvAp_b);

    } else {    // ~correctA
        if (scale) {
            auto sEpvAp(diag_matrix(wEpvAp_b,'-'));

            ret->M.reset(new EigenSparseMatrix(sEpvAp * *EpvAp));
        } else {
            ret->M = std::move(EpvAp);
        }

        // ----- Compute the final weight (diagonal) matrix
        ret->weight.reference(wEpvAp_b);
    }

    return ret;
}
// ----------------------------------------------------------------
RegridMatrices::RegridMatrices(IceRegridder *regridder)
{
    printf("===== RegridMatrices Grid geometries:\n");
    printf("    nA = %d\n", regridder->gcm->nA());
    printf("    nhc = %d\n", regridder->gcm->nhc());
    printf("    nE = %d\n", regridder->gcm->nE());
    printf("    nI = %d\n", regridder->nI());
    printf("    nG = %d\n", regridder->nG());

    UrAE urA(regridder->gcm->nA(),
        std::bind(&IceRegridder::GvAp, regridder, _1),
        std::bind(&IceRegridder::sApvA, regridder, _1));

    UrAE urE(regridder->gcm->nE(),
        std::bind(&IceRegridder::GvEp, regridder, _1),
        std::bind(&IceRegridder::sEpvE, regridder, _1));

    // ------- AvI, IvA
    regrids.insert(make_pair("AvI", std::bind(&compute_AEvI, regridder, _1, _2, urA) ));
    regrids.insert(make_pair("IvA", std::bind(&compute_IvAE, regridder, _1, _2, urA) ));

    // ------- EvI, IvE
    regrids.insert(make_pair("EvI", std::bind(&compute_AEvI, regridder, _1, _2, urE) ));
    regrids.insert(make_pair("IvE", std::bind(&compute_IvAE, regridder, _1, _2, urE) ));

    // ------- EvA, AvE regrids.insert(make_pair("EvA", std::bind(&compute_EvA, this, _1, _2, urE, urA) ));
    regrids.insert(make_pair("EvA", std::bind(&compute_EvA, regridder, _1, _2, urE, urA) ));
    regrids.insert(make_pair("AvE", std::bind(&compute_EvA, regridder, _1, _2, urA, urE) ));

    // ----- Show what we have!
    printf("Available Regrids:");
    for (auto ii = regrids.begin(); ii != regrids.end(); ++ii) {
        printf(" %s", ii->first.c_str());
    }
    printf("\n");
}


} // namespace
