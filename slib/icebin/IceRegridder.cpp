#include <cstdio>
#include <iostream>
#include <icebin/GCMRegridder.hpp>
#include <icebin/IceRegridder_L0.hpp>
#include <icebin/smoother.hpp>
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
void IceRegridder::set_elevI(DenseArrayT<1> const &_elevI)
{
    elevI = _elevI;        // Copies
}
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
void IceRegridder::sApvA(MakeDenseEigenT::AccumT &w)
{
    ibmisc::Proj_LL2XY proj(gridI->sproj);
    for (auto cell=gcm->gridA->cells.begin(); cell != gcm->gridA->cells.end(); ++cell) {
        w.add({cell->index, cell->index}, cell->native_area / cell->proj_area(&proj));
    }
}

/** Produces the diagonal matrix [Elevation projected] <-- [Elevation]
NOTE: wAvAp == sApvA */
void IceRegridder::sEpvE(MakeDenseEigenT::AccumT &w)
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
    if (rw_full) ncio_blitz(ncio, elevI, true, vname + ".elevI", "double",
        get_dims(ncio ,{vname + ".gridI.cells.nfull"}));

printf("END IceRegridder::ncio(%s, %d)\n", vname.c_str(), rw_full);
}

void IceRegridder::init(
    std::string const &name,
    std::unique_ptr<Grid> &&_gridI,
    std::unique_ptr<Grid> &&_exgrid,
    InterpStyle _interp_style,
    blitz::Array<double,1> const &_elevI)
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

    typedef std::function<void(MakeDenseEigenT::AccumT &)> matrix_fn;

    const matrix_fn GvAp;
    const matrix_fn sApvA;

    UrAE(long _nfull, matrix_fn _GvAp, matrix_fn _sApvA) : nfull(_nfull), GvAp(_GvAp), sApvA(_sApvA) {}
};



// ------------------------------------------------------------
std::function<bool(long)> in_sparse_fn(SparseSetT const &dim)
    { return std::bind(&SparseSetT::in_sparse, &dim, _1); }
// ------------------------------------------------------------
static std::unique_ptr<WeightedSparse> compute_AEvI(IceRegridder *regridder,
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &params, UrAE const &AE)
{
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse(dims));
    SparseSetT * const dimA(ret->dims[0]);
    SparseSetT * const dimI(ret->dims[1]);
    SparseSetT _dimG;
    SparseSetT * const dimG(&_dimG);

    if (dimA) dimA->set_sparse_extent(AE.nfull);
    if (dimI) dimI->set_sparse_extent(regridder->nI());
    dimG->set_sparse_extent(regridder->nG());

    // ----- Get the Ur matrices (which determines our dense dimensions)
    MakeDenseEigenT GvI_m(
        std::bind(&IceRegridder::GvI, regridder, _1),
        SparsifyTransform::ADD_DENSE,
        {dimG, dimI}, '.');
    MakeDenseEigenT ApvG_m(        // _m ==> type MakeDenseEigenT
        AE.GvAp,
        SparsifyTransform::ADD_DENSE,
        {dimG, dimA}, 'T');

    // ----- Convert to Eigen and multiply
    auto ApvG(ApvG_m.to_eigen());
    auto GvI(GvI_m.to_eigen());
    auto sGvI(sum_to_diagonal(GvI, 0, '-'));

    std::unique_ptr<EigenSparseMatrixT> ApvI(
        new EigenSparseMatrixT(ApvG * sGvI * GvI));
    ret->Mw.reference(sum(*ApvI, 1, '+'));    // Area of I cells

    // ----- Apply final scaling, and convert back to sparse dimension
    if (params.correctA) {
        // ----- Compute the final weight matrix
        auto wAvAp(MakeDenseEigenT(                   // diagonal
            AE.sApvA,
            SparsifyTransform::TO_DENSE_IGNORE_MISSING,
            {dimA, dimA}, '.').to_eigen());
        auto wApvI(sum_to_diagonal(*ApvI, 0, '+'));        // diagonal
        EigenSparseMatrixT wAvI(wAvAp * wApvI);    // diagonal...

        // +correctA: Weight matrix in A space
        ret->wM.reference(sum(wAvI, 0, '+'));    // Area of A cells

        // Compute the main matrix
        if (params.scale) {
            // Get two diagonal Eigen scale matrices
            auto sAvAp(sum_to_diagonal(wAvAp, 0, '-'));
            auto sApvI(sum_to_diagonal(wApvI, 0, '-'));

            ret->M.reset(new EigenSparseMatrixT(
                sAvAp * sApvI * *ApvI));    // AvI_scaled
        } else {
            ret->M = std::move(ApvI);
        }

    } else {

        // ----- Compute the final weight matrix
        // ~correctA: Weight matrix in Ap space
        auto wApvI_b(sum(*ApvI, 0, '+'));
        ret->wM.reference(wApvI_b);

        if (params.scale) {
            // Get two diagonal Eigen scale matrices
            auto sApvI(diag_matrix(wApvI_b, '-'));

            ret->M.reset(new EigenSparseMatrixT(
                sApvI * *ApvI));    // ApvI_scaled
        } else {
            ret->M = std::move(ApvI);
        }
    }

    return ret;
}
// ---------------------------------------------------------
// ---------------------------------------------------------
std::unique_ptr<WeightedSparse> compute_IvAE(IceRegridder *regridder,
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &params, UrAE const &AE)
{
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse(dims));
    SparseSetT * const dimA(ret->dims[1]);
    SparseSetT * const dimI(ret->dims[0]);
    SparseSetT _dimG;
    SparseSetT * const dimG(&_dimG);

    if (dimA) dimA->set_sparse_extent(AE.nfull);
    if (dimI) dimI->set_sparse_extent(regridder->nI());
    dimG->set_sparse_extent(regridder->nG());

    // ----- Get the Ur matrices (which determines our dense dimensions)
    MakeDenseEigenT GvAp_m(
        AE.GvAp,
        SparsifyTransform::ADD_DENSE,
        {dimG, dimA}, '.');
    MakeDenseEigenT IvG_m(
        std::bind(&IceRegridder::GvI, regridder, _1),
        SparsifyTransform::ADD_DENSE,
        {dimG, dimI}, 'T');

    // ----- Convert to Eigen
    auto GvAp(GvAp_m.to_eigen());
    auto IvG(IvG_m.to_eigen());
    auto sGvAp(sum_to_diagonal(GvAp, 0, '-'));

    // Unscaled matrix
    std::unique_ptr<EigenSparseMatrixT> IvAp(
        new EigenSparseMatrixT(IvG * sGvAp * GvAp));

    // Get weight vector from IvAp_e
    ret->wM.reference(sum(*IvAp, 0, '+'));

    // ----- Apply final scaling, and convert back to sparse dimension
    if (params.correctA) {
        // Scaling matrix
        auto sApvA(MakeDenseEigenT(
            AE.sApvA,
            SparsifyTransform::TO_DENSE_IGNORE_MISSING,
            {dimA, dimA}, '.').to_eigen());

        // Compute area of A grid cells
        auto IvApw(sum_to_diagonal(*IvAp, 1, '+'));    // Area of A cells
        auto &wAvAp(sApvA);    // Symmetry: wAvAp == sApvA
        EigenSparseMatrixT Aw(wAvAp * IvApw);
        ret->Mw.reference(sum(Aw,0,'+'));

        if (params.scale) {
            auto sIvAp(sum_to_diagonal(*IvAp, 0, '-'));
            ret->M.reset(new EigenSparseMatrixT(
                sIvAp * *IvAp * sApvA));
        } else {
            ret->M.reset(new EigenSparseMatrixT(
                *IvAp * sApvA));
        }
    } else {
        ret->Mw.reference(sum(*IvAp, 1, '+'));    // Area of A cells
        if (params.scale) {
            auto sIvAp(sum_to_diagonal(*IvAp, 0, '-'));
            ret->M.reset(new EigenSparseMatrixT(
                sIvAp * *IvAp));
        } else {
            ret->M = std::move(IvAp);
        }
    }

    return ret;
}

static std::unique_ptr<WeightedSparse> compute_EvA(IceRegridder *regridder,
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &params, UrAE const &E, UrAE const &A)
{
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse(dims));
    SparseSetT * const dimE(ret->dims[0]);
    SparseSetT * const dimA(ret->dims[1]);
    SparseSetT _dimG;
    SparseSetT * const dimG(&_dimG);

    if (dimA) dimA->set_sparse_extent(A.nfull);
    if (dimE) dimE->set_sparse_extent(E.nfull);
    dimG->set_sparse_extent(regridder->nG());

    // ----- Get the Ur matrices (which determines our dense dimensions)

    MakeDenseEigenT GvAp_m(
        A.GvAp,
        SparsifyTransform::ADD_DENSE,
        {dimG, dimA}, '.');
    MakeDenseEigenT EpvG_m(
        E.GvAp,
        SparsifyTransform::ADD_DENSE,
        {dimG, dimE}, 'T');

    // ----- Convert to Eigen and multiply
    auto GvAp(GvAp_m.to_eigen());
    auto EpvG(EpvG_m.to_eigen());

    auto sGvAp(sum_to_diagonal(GvAp, 0, '-'));

    // Unweighted matrix
    std::unique_ptr<EigenSparseMatrixT> EpvAp(
        new EigenSparseMatrixT(EpvG * sGvAp * GvAp));

    // ----- Apply final scaling, and convert back to sparse dimension
    auto wEpvAp(sum_to_diagonal(*EpvAp,0,'+'));
    if (params.correctA) {
        auto sApvA(MakeDenseEigenT(
            A.sApvA,
            SparsifyTransform::TO_DENSE_IGNORE_MISSING,
            {dimA, dimA}, '.').to_eigen());

        auto wEvEp(MakeDenseEigenT(
            E.sApvA,
            SparsifyTransform::TO_DENSE_IGNORE_MISSING,
            {dimE, dimE}, '.').to_eigen());

        // +correctA: Weight matrix in E space
        EigenSparseMatrixT wEvAp(wEvEp * wEpvAp);
        ret->wM.reference(sum(wEvAp,0,'+'));

        // Compute area of A cells
        auto EpvApw(sum_to_diagonal(*EpvAp,1,'+'));
        auto &wAvAp(sApvA);    // Symmetry: wAvAp == sApvA
        EigenSparseMatrixT Aw(wAvAp * EpvApw);
        ret->Mw.reference(sum(Aw,0,'+'));

        if (params.scale) {
            auto sEvAp(sum_to_diagonal(wEvAp,0,'-'));
            ret->M.reset(new EigenSparseMatrixT(
                sEvAp * *EpvAp * sApvA));    // EvA
        } else {
            ret->M.reset(new EigenSparseMatrixT(*EpvAp * sApvA));
        }
    } else {    // ~correctA
        // ~correctA: Weight matrix in Ep space
        ret->wM.reference(sum(wEpvAp,0,'+'));
        ret->Mw.reference(sum(*EpvAp,1,'+'));
        if (params.scale) {
            auto sEpvAp(sum_to_diagonal(wEpvAp,0,'-'));

            ret->M.reset(new EigenSparseMatrixT(sEpvAp * *EpvAp));
        } else {
            ret->M = std::move(EpvAp);
        }
    }

    return ret;
}
// ----------------------------------------------------------------
void RegridMatrices::add_regrid(std::string const &spec,
    RegridMatrices::SmoothingFunction const *smooth,
    RegridMatrices::MatrixFunction const &regrid)
{
    regrids.insert(make_pair(spec, Binding(regrid, smooth)));
}
// ----------------------------------------------------------------
RegridMatrices::RegridMatrices(IceRegridder *regridder)
{
#if 0
    printf("===== RegridMatrices Grid geometries:\n");
    printf("    nA = %d\n", regridder->gcm->nA());
    printf("    nhc = %d\n", regridder->gcm->nhc());
    printf("    nE = %d\n", regridder->gcm->nE());
    printf("    nI = %d\n", regridder->nI());
    printf("    nG = %d\n", regridder->nG());
#endif

    UrAE urA(regridder->gcm->nA(),
        std::bind(&IceRegridder::GvAp, regridder, _1),
        std::bind(&IceRegridder::sApvA, regridder, _1));

    UrAE urE(regridder->gcm->nE(),
        std::bind(&IceRegridder::GvEp, regridder, _1),
        std::bind(&IceRegridder::sEpvE, regridder, _1));

    // Get a smoothing function for each output grid
     smoothI = std::bind(&smoothing_matrix,
        _1, &*regridder->gridI, _2, &regridder->elevI, _3, _4);

    // ------- AvI, IvA
    add_regrid("AvI", nullptr,
        std::bind(&compute_AEvI, regridder, _1, _2, urA));
    add_regrid("IvA", &smoothI,
        std::bind(&compute_IvAE, regridder, _1, _2, urA));

    // ------- EvI, IvE
    add_regrid("EvI", nullptr,
        std::bind(&compute_AEvI, regridder, _1, _2, urE));
    add_regrid("IvE", &smoothI,
        std::bind(&compute_IvAE, regridder, _1, _2, urE));

    // ------- EvA, AvE regrids.insert(make_pair("EvA", std::bind(&compute_EvA, this, _1, _2, urE, urA) ));
    add_regrid("EvA", nullptr,
        std::bind(&compute_EvA, regridder, _1, _2, urE, urA));
    add_regrid("AvE", nullptr,
        std::bind(&compute_EvA, regridder, _1, _2, urA, urE));

#if 0
    // ----- Show what we have!
    printf("Available Regrids:");
    for (auto ii = regrids.begin(); ii != regrids.end(); ++ii) {
        printf(" %s", ii->first.c_str());
    }
    printf("\n");
#endif

}
// ----------------------------------------------------------------
std::unique_ptr<WeightedSparse> RegridMatrices::matrix(
    std::string const &spec_name,
    std::array<SparseSetT *,2> dims,
    Params const &params) const
{
    auto &bindings(regrids.at(spec_name));
    auto BvA(bindings.matrix(dims, params));
    // If the user asked for smoothing and we can smooth this kind of matrix...
    if (params.smooth() && bindings.smooth) {
        TupleListT<2> smoothB_t({dims[0]->dense_extent(), dims[0]->dense_extent()});
        (*bindings.smooth)(smoothB_t, *dims[0], BvA->wM, params.sigma);
        auto smoothB(to_eigen(smoothB_t));

        // Smooth the underlying unsmoothed regridding transformation
        BvA->M.reset(new EigenSparseMatrixT(smoothB * *BvA->M));
        BvA->smooth = true;
        BvA->conserve = params.conserve;
    }
    return BvA;
}
// ----------------------------------------------------------------
static void mask_result(EigenDenseMatrixT &ret, blitz::Array<double,1> const &wB_b, double fill)
{
    int nB = ret.rows();    // == wB_b.extent(0)
    int nvar = ret.cols();

    // Mask out cells that slipped into the output because they were
    // in the SparseSet; but don't actually get any contribution.
    for (int i=0; i<nB; ++i) {
        if (wB_b(i) != 0.) continue;
        for (int n=0; n<nvar; ++n) ret(i,n) = fill;
    }

}

EigenDenseMatrixT WeightedSparse::apply_e(
    // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
    blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable
    double fill) const    // Fill value for cells not in BvA matrix
{
    auto &BvA(*this);

    // A{jn}   One col per variable
    int nvar = A_b.extent(0);
    int nA = A_b.extent(1);
    Eigen::Map<EigenDenseMatrixT> const A(
        const_cast<double *>(A_b.data()), nA, nvar);

    // |i| = size of output vector space (B)
    // |j| = size of input vector space (A)
    // |n| = number of variables being processed together

    // Apply initial regridding.
    EigenDenseMatrixT B0(*BvA.M * A);        // B0{in}

    // Only apply conservation correction if all of:
    //   a) Matrix is smoothed, so it needs a conservation correction
    //   b) User requested conservation be maintained
    if (!(BvA.smooth && BvA.conserve)) {
        // Remove cells not in the sparse matrix
        mask_result(B0, BvA.wM, fill);
        return B0;
    }

    // -------------- Apply the Conservation Correction
    // Integrate each variable of input (A) over full domain
    auto &wA_b(BvA.Mw);
    Eigen::Map<EigenRowVectorT> const wA(const_cast<double *>(wA_b.data()), 1, wA_b.extent(0));
    typedef Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic> EigenArrayT;
    EigenArrayT TA((wA * A).array());        // TA{n} row array

    // Integrate each variable of output (B) over full domain
    auto &wB_b(BvA.wM);
    int nB = wB_b.extent(0);
    Eigen::Map<EigenRowVectorT> const wB(const_cast<double *>(wB_b.data()), 1, wB_b.extent(0));
    EigenArrayT TB((wB * B0).array());
    EigenArrayT TB_inv(1. / TB);    // TB_inv{n} row array

    // Factor{nn}: Conservation correction for each variable.

    auto Factor(TA * TB_inv);    // Factor{n} row array

    std::cout << "-------- WeightedSparse::apply() conservation" << std::endl;
    std::cout << "    |input|    = " << TA << std::endl;
    std::cout << "    |output|   = " << TB << std::endl;
    std::cout << "    correction = " << Factor << std::endl;

    EigenDenseMatrixT ret(B0 * Factor.matrix().asDiagonal());    // ret{in}
    // Remove cells not in the sparse matrix
    mask_result(B0, BvA.wM, fill);

    return ret;
}

void WeightedSparse::ncio(ibmisc::NcIO &ncio,
    std::string const &vname,
    std::array<std::string,2> dim_names)
{

    auto ncdims(ibmisc::get_or_add_dims(ncio,
        {dim_names[0] + ".dense_extent", dim_names[1] + ".dense_extent"},
        {dims[0]->dense_extent(), dims[1]->dense_extent()}));

    // ----------- wM
    std::string matrix_name(dim_names[0] + "v" + dim_names[1]);
    ncio_blitz<double,1>(ncio, wM, true,
        vname + ".wM",
        {ncdims[0]});

    // --------- M
#if 1
    ncio_eigen(ncio, *M,
        vname + ".M");

    netCDF::NcVar ncvar = ncio.nc->getVar(vname + ".M.info");
    get_or_put_att(ncvar, ncio.rw,
        "smooth", get_nc_type<bool>(), &smooth, 1);
    get_or_put_att(ncvar, ncio.rw,
        "conserve", get_nc_type<bool>(), &conserve, 1);
#endif

    // ---- Mw
    ncio_blitz<double,1>(ncio, Mw, true,
        vname + ".Mw",
        {ncdims[1]});
}


} // namespace
