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
void IceRegridder::sApvA(SparseTriplets<SparseMatrix> &w, std::function<bool(long)> const &filter_fn)
{
    ibmisc::Proj_LL2XY proj(gridI->sproj);
    for (auto cell=gcm->gridA->cells.begin(); cell != gcm->gridA->cells.end(); ++cell) {
        if (!filter_fn(cell->index)) continue;

        w.add({cell->index, cell->index}, cell->native_area / cell->proj_area(&proj));
    }
}

/** Produces the diagonal matrix [Elevation projected] <-- [Elevation]
NOTE: wAvAp == sApvA */
void IceRegridder::sEpvE(SparseTriplets<SparseMatrix> &w, std::function<bool(long)> const &filter_fn)
{
    ibmisc::Proj_LL2XY proj(gridI->sproj);
    for (auto cell=gcm->gridA->cells.begin(); cell != gcm->gridA->cells.end(); ++cell) {
        long nhp = gcm->nhp(cell->index);
        long tuple[2] = {cell->index, 0};
        long &ihp(tuple[1]);
        for (ihp=0; ihp<nhp; ++ihp) {
            long indexE = gcm->indexingHC.tuple_to_index(tuple);
            if (!filter_fn(indexE)) continue;
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
    weights[0] = (1.0 - ratio);
    weights[1] = ratio;
}
// ========================================================================

/** Used to generate Ur matrices for Atm or Elevation grids.
This allows us to exploit an algebraic symmetry between A and E. */
struct UrAE {
    const long nfull;

    typedef std::function<void(spsparse::SparseTriplets<SparseMatrix> &)> GvAp_fn;
    const GvAp_fn GvAp;

    typedef const std::function<void(
        spsparse::SparseTriplets<SparseMatrix> &w,
        std::function<bool(long)>
    )> sApvA_fn;
    const sApvA_fn sApvA;

    UrAE(long _nfull, GvAp_fn _GvAp, sApvA_fn _sApvA) : nfull(_nfull), GvAp(_GvAp), sApvA(_sApvA) {}
};



// ------------------------------------------------------------
std::function<bool(long)> in_sparse_fn(SparseSet const &dim)
    { return std::bind(&SparseSet::in_sparse, &dim, _1); }
// ------------------------------------------------------------
static std::unique_ptr<WeightedSparse> compute_AEvI(IceRegridder *regridder, bool scale, bool correctA, UrAE const &AE)
{
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse);
    auto &dimA(ret->dims[0]);
    auto &dimI(ret->dims[1]);
    SparseSet dimG;

    // ----- Get the Ur matrices (which determines our dense dimensions)
    // GvI
    SparseTriplets<SparseMatrix> GvI_t({&dimG, &dimI});    // _t=triplets
    GvI_t.set_shape({regridder->nG(), regridder->nI()});
    regridder->GvI(GvI_t);

    // ApvG
    SparseTriplets<SparseMatrix> GvAp_t({&dimG, &dimA});
    GvAp_t.set_shape({regridder->nG(), AE.nfull});
    AE.GvAp(GvAp_t);

    // ----- Convert to Eigen and multiply
    auto GvI_e(GvI_t.to_eigen());
    auto sGvI_e(scale_matrix(GvI_e,0));
    auto ApvG_e(GvAp_t.to_eigen('T'));
    std::unique_ptr<EigenSparseMatrix> ApvI_e(
        new EigenSparseMatrix(ApvG_e * sGvI_e * GvI_e));

    // ----- Apply final scaling, and convert back to sparse dimension
    if (correctA) {
        // Scaling matrix AvAp
        SparseTriplets<SparseMatrix> sApvA({&dimA, &dimA});
        sApvA.set_shape({AE.nfull, AE.nfull});
        AE.sApvA(sApvA, in_sparse_fn(dimA));   // Filter to avoid adding more items to dimA

        // ----- Compute the final weight matrix
        auto wAvAp_e(sApvA.to_eigen('.', false));       // Invert it...
        auto wApvI_e(weight_matrix(*ApvI_e, 0));
        EigenSparseMatrix wAvI_e(wAvAp_e * wApvI_e);    // diagonal...

        // Sums rows of an Eigen matrix into a dense blitz::Array
        ret->weight.reference(sum(wAvI_e,0));

        // Compute the main matrix
        if (scale) {
            // Get two diagonal Eigen scale matrices
            auto sAvAp_e(sApvA.to_eigen('.', true));        // Invert it...
            auto sApvI_e(scale_matrix(*ApvI_e, 0));

            ret->M.reset(new EigenSparseMatrix(
                sAvAp_e * sApvI_e * *ApvI_e));    // AvI_scaled
        } else {
            ret->M = std::move(ApvI_e);
        }

    } else {

        // ----- Compute the final weight matrix
        ret->weight.reference(sum(*ApvI_e, 0));    // wApvI

        if (scale) {
            // Get two diagonal Eigen scale matrices
            auto sApvI_e(scale_matrix(*ApvI_e, 0));

            ret->M.reset(new EigenSparseMatrix(
                sApvI_e * *ApvI_e));    // ApvI_scaled
        } else {
            ret->M = std::move(ApvI_e);
        }

    }

    return ret;
}

static std::unique_ptr<WeightedSparse> compute_IvAE(IceRegridder *regridder, bool scale, bool correctA, UrAE const &AE)
{
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse);
    auto &dimA(ret->dims[1]);
    auto &dimI(ret->dims[0]);
    SparseSet dimG;

    // ----- Get the Ur matrices (which determines our dense dimensions)
    // ApvG
    SparseTriplets<SparseMatrix> GvAp({&dimG, &dimA});
    GvAp.set_shape({regridder->nG(), AE.nfull});
    AE.GvAp(GvAp);

    // GvI
    SparseTriplets<SparseMatrix> GvI({&dimG, &dimI});
    GvI.set_shape({regridder->nG(), regridder->nI()});
    regridder->GvI(GvI);

    // ----- Convert to Eigen and multiply
    auto GvAp_e(GvAp.to_eigen('.'));
    auto sGvAp_e(GvAp.eigen_scale_matrix(0));
    auto IvG_e(GvI.to_eigen('T'));

    // Unscaled matrix
    std::unique_ptr<EigenSparseMatrix> IvAp_e(
        new EigenSparseMatrix(IvG_e * sGvAp_e * GvAp_e));

    // Get weight vector from IvAp_e
    ret->weight.reference(sum(*IvAp_e, 0));

    // ----- Apply final scaling, and convert back to sparse dimension
    if (correctA) {
        // sApvA
        SparseTriplets<SparseMatrix> sApvA({&dimA, &dimA});
        sApvA.set_shape({AE.nfull, AE.nfull});
        AE.sApvA(sApvA, in_sparse_fn(dimA)); // Filter to avoid adding more items to dimA
        auto sApvA_e(sApvA.to_eigen('.', false));       // Don't invert it...

        if (scale) {
            auto sIvAp_e(scale_matrix(*IvAp_e, 0));
            ret->M.reset(new EigenSparseMatrix(sIvAp_e * *IvAp_e * sApvA_e));
        } else {
            ret->M.reset(new EigenSparseMatrix(*IvAp_e * sApvA_e));
        }
    } else {
        if (scale) {
            auto sIvAp_e(scale_matrix(*IvAp_e, 0));
            ret->M.reset(new EigenSparseMatrix(sIvAp_e * *IvAp_e));
        } else {
            ret->M = std::move(IvAp_e);
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

    // ApvG
    SparseTriplets<SparseMatrix> GvAp({&dimG, &dimA});
    GvAp.set_shape({regridder->nG(), A.nfull});
    A.GvAp(GvAp);

    // GvEp
    SparseTriplets<SparseMatrix> GvEp({&dimG, &dimE});
    GvEp.set_shape({regridder->nG(), E.nfull});
    E.GvAp(GvEp);


    // ----- Convert to Eigen and multiply
    auto GvAp_e(GvAp.to_eigen('.'));
    auto sGvAp_e(GvAp.eigen_scale_matrix(0));
    auto EpvG_e(GvEp.to_eigen('T'));

    // Unweighted matrix
    std::unique_ptr<EigenSparseMatrix> EpvAp_e(
        new EigenSparseMatrix(EpvG_e * sGvAp_e * GvAp_e));

    // ----- Apply final scaling, and convert back to sparse dimension
    if (correctA) {

        // sApvA
        SparseTriplets<SparseMatrix> sApvA_t({&dimA, &dimA});
        sApvA_t.set_shape({A.nfull, A.nfull});
        A.sApvA(sApvA_t, in_sparse_fn(dimA));  // Filter to avoid adding more items to dimA
        auto sApvA_e(sApvA_t.to_eigen('.', false));       // Don't invert it...

        // ---- Obtain sEpvE
        SparseTriplets<SparseMatrix> sEpvE_t({&dimE, &dimE});
        sEpvE_t.set_shape({E.nfull, E.nfull});
        E.sApvA(sEpvE_t, std::bind(&SparseSet::in_sparse, &dimE, _1));    // Avoid adding more items to dimE

        if (scale) {
            auto sEvEp_e(sEpvE_t.to_eigen('.', true));        // Invert it...
            auto sEpvAp_e(scale_matrix(*EpvAp_e, 0));

            ret->M.reset(new EigenSparseMatrix(
                sEvEp_e * sEpvAp_e * *EpvAp_e * sApvA_e));    // EvA
        } else {
            ret->M.reset(new EigenSparseMatrix(*EpvAp_e * sApvA_e));
        }

        // ----- Compute the final weight (diagonal) matrix
        auto wEvEp_e(sEpvE_t.to_eigen('.', false));
        auto wEpvAp_b(sum(*EpvAp_e,0));    // _b = Blitz++
        EigenSparseMatrix weight_e(wEvEp_e * weight_matrix(wEpvAp_b));
        ret->weight.reference(wEpvAp_b);

    } else {    // ~correctA
        auto wEpvAp_b(sum(*EpvAp_e,0));
        if (scale) {
            auto sEpvAp_e(scale_matrix(wEpvAp_b));

            ret->M.reset(new EigenSparseMatrix(sEpvAp_e * *EpvAp_e));
        } else {
            ret->M = std::move(EpvAp_e);
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
    printf("    nhp = %d\n", regridder->gcm->nhp());
    printf("    nE = %d\n", regridder->gcm->nE());
    printf("    nI = %d\n", regridder->nI());
    printf("    nG = %d\n", regridder->nG());

    UrAE urA(regridder->gcm->nA(),
        std::bind(&IceRegridder::GvAp, regridder, _1),
        std::bind(&IceRegridder::sApvA, regridder, _1, _2));

    UrAE urE(regridder->gcm->nE(),
        std::bind(&IceRegridder::GvEp, regridder, _1),
        std::bind(&IceRegridder::sEpvE, regridder, _1, _2));

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
