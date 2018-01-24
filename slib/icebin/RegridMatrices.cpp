#ifndef ICEBIN_REGRID_MATRICES_CPP
#define ICEBIN_REGRID_MATRICES_CPP

#include <functional>

#include <icebin/RegridMatrices.hpp>
#include <icebin/IceRegridder.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/smoother.hpp>
#include <icebin/IceRegridder_L0.hpp>

using namespace std::placeholders;
using namespace spsparse;
using namespace ibmisc;

namespace icebin {

// =======================================================================
// Regrid Matrix Generation

/** Used to generate Ur matrices for Atm or Elevation grids.
This allows us to exploit an algebraic symmetry between A and E. */
struct UrAE {
    std::string name;    // For deugging
    const long nfull;

    typedef std::function<void(MakeDenseEigenT::AccumT &&)> matrix_fn;

    const matrix_fn GvAp;
    const matrix_fn sApvA;

    UrAE(std::string const &_name, long _nfull, matrix_fn _GvAp, matrix_fn _sApvA) : name(_name), nfull(_nfull), GvAp(_GvAp), sApvA(_sApvA) {}
};



// ------------------------------------------------------------
static std::unique_ptr<WeightedSparse> compute_AEvI(
    IceRegridder const *regridder,
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &params,
    blitz::Array<double,1> const *elevI,
    UrAE const &AE)
{
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse(dims, true));
    SparseSetT * const dimA(ret->dims[0]);
    SparseSetT * const dimI(ret->dims[1]);
    SparseSetT _dimG;
    SparseSetT * const dimG(&_dimG);

    if (dimA) dimA->set_sparse_extent(AE.nfull);
    if (dimI) dimI->set_sparse_extent(regridder->nI());
    dimG->set_sparse_extent(regridder->nG());

    // ----- Get the Ur matrices (which determines our dense dimensions)
    MakeDenseEigenT GvI_m(
        // Only includes ice model grid cells with ice in them.
        std::bind(&IceRegridder::GvI, regridder, _1, elevI),
        {SparsifyTransform::ADD_DENSE},
        {dimG, dimI}, '.');
    MakeDenseEigenT ApvG_m(        // _m ==> type MakeDenseEigenT
        // Only includes ice model grid cells with ice in them.
        AE.GvAp,
        {SparsifyTransform::ADD_DENSE},
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
            {SparsifyTransform::TO_DENSE_IGNORE_MISSING},
            {dimA, dimA}, '.').to_eigen());
        auto wApvI(sum_to_diagonal(*ApvI, 0, '+'));        // diagonal
        EigenSparseMatrixT wAvI(wAvAp * wApvI);    // diagonal...

        // +correctA: Weight matrix in A space
        ret->wM.reference(sum(wAvI, 0, '+'));    // Area of A cells

        // Compute the main matrix
        auto sAvAp(sum_to_diagonal(wAvAp, 0, '-'));
        if (params.scale) {
            // Get two diagonal Eigen scale matrices
            auto sApvI(sum_to_diagonal(wApvI, 0, '-'));

            ret->M.reset(new EigenSparseMatrixT(
                sAvAp * sApvI * *ApvI));    // AvI_scaled
        } else {
            // Should be like this for test_conserv.py
            // Note that sAvAp * sApvI = [size (weight) of grid cells in A]
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
std::unique_ptr<WeightedSparse> compute_IvAE(
    IceRegridder const *regridder,
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &params,
    blitz::Array<double,1> const *elevI,
    UrAE const &AE)
{
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse(dims, !params.smooth()));
    SparseSetT * const dimA(ret->dims[1]);    SparseSetT * const dimI(ret->dims[0]);
    SparseSetT _dimG;
    SparseSetT * const dimG(&_dimG);

    if (dimA) dimA->set_sparse_extent(AE.nfull);
    if (dimI) dimI->set_sparse_extent(regridder->nI());
    dimG->set_sparse_extent(regridder->nG());

    // ----- Get the Ur matrices (which determines our dense dimensions)
    MakeDenseEigenT GvAp_m(
        AE.GvAp,
        {SparsifyTransform::ADD_DENSE},
        {dimG, dimA}, '.');
    MakeDenseEigenT IvG_m(
        std::bind(&IceRegridder::GvI, regridder, _1, elevI),
        {SparsifyTransform::ADD_DENSE},
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
            {SparsifyTransform::TO_DENSE_IGNORE_MISSING},
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

    // Smooth the result on I, if needed
    if (params.smooth()) {

        // Obtain the smoothing matrix (smoother.hpp)
        TupleListT<2> smoothI_t({dimI->dense_extent(), dimI->dense_extent()});
        smoothing_matrix(smoothI_t, regridder->agridI,
            *dimI, *elevI, ret->wM, params.sigma);
        auto smoothI(to_eigen_sparsematrix(smoothI_t));

        // Smooth the underlying unsmoothed regridding transformation
        ret->M.reset(new EigenSparseMatrixT(smoothI * *ret->M));
    }
    return ret;
}

static std::unique_ptr<WeightedSparse> compute_EvA(IceRegridder const *regridder,
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &params, UrAE const &E, UrAE const &A)
{
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse(dims, true));
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
        {SparsifyTransform::ADD_DENSE},
        {dimG, dimA}, '.');
    MakeDenseEigenT EpvG_m(
        E.GvAp,
        {SparsifyTransform::ADD_DENSE},
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
            {SparsifyTransform::TO_DENSE_IGNORE_MISSING},
            {dimA, dimA}, '.').to_eigen());

        auto wEvEp(MakeDenseEigenT(
            E.sApvA,
            {SparsifyTransform::TO_DENSE_IGNORE_MISSING},
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


RegridMatrices GCMRegridder_Standard::regrid_matrices(
    int sheet_index,
    blitz::Array<double,1> const &elevI) const
{
    IceRegridder const *regridder = &*ice_regridders()[sheet_index];

#if 0
    printf("===== RegridMatrices Grid geometries:\n");
    printf("    nA = %d\n", this->nA());
    printf("    nhc = %d\n", this->nhc());
    printf("    nE = %d\n", this->nE());
    printf("    nI = %d\n", regridder->nI());
    printf("    nG = %d\n", regridder->nG());
#endif

    RegridMatrices rm(regridder);

    UrAE urA("UrA", this->nA(),
        std::bind(&IceRegridder::GvAp, regridder, _1, &elevI),
        std::bind(&IceRegridder::sApvA, regridder, _1));

    UrAE urE("UrE", this->nE(),
        std::bind(&IceRegridder::GvEp, regridder, _1, &elevI),
        std::bind(&IceRegridder::sEpvE, regridder, _1));

    // ------- AvI, IvA
    rm.add_regrid("AvI",
        std::bind(&compute_AEvI, regridder, _1, _2, &elevI, urA));
    rm.add_regrid("IvA",
        std::bind(&compute_IvAE, regridder, _1, _2, &elevI, urA));

    // ------- EvI, IvE
    rm.add_regrid("EvI",
        std::bind(&compute_AEvI, regridder, _1, _2, &elevI, urE));
    rm.add_regrid("IvE",
        std::bind(&compute_IvAE, regridder, _1, _2, &elevI, urE));

    // ------- EvA, AvE regrids.insert(make_pair("EvA", std::bind(&compute_EvA, regridder, _1, _2, urE, urA) ));
    rm.add_regrid("EvA",
        std::bind(&compute_EvA, regridder, _1, _2, urE, urA));
    rm.add_regrid("AvE",
        std::bind(&compute_EvA, regridder, _1, _2, urA, urE));

#if 0
    // ----- Show what we have!
    printf("Available Regrids:");
    for (auto ii = regrids.begin(); ii != regrids.end(); ++ii) {
        printf(" %s", ii->first.c_str());
    }
    printf("\n");
#endif

    return std::move(rm);
}
// -----------------------------------------------------------------------
// ----------------------------------------------------------------
void RegridMatrices::add_regrid(std::string const &spec,
    RegridMatrices::MatrixFunction const &regrid)
{
    regrids.insert(make_pair(spec, regrid));
}
// ----------------------------------------------------------------
// ----------------------------------------------------------------
std::unique_ptr<WeightedSparse> RegridMatrices::matrix(
    std::string const &spec_name,
    std::array<SparseSetT *,2> dims,
    Params const &params) const
{
    auto &matrix_fn(regrids.at(spec_name));
    auto BvA(matrix_fn(dims, params));
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
// -----------------------------------------------------------------------
EigenDenseMatrixT WeightedSparse::apply_e(
    // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
    blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable
    double fill,     // Fill value for cells not in BvA matrix
    bool force_conservation) const
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


    if (BvA.M->cols() != A.rows()) (*icebin_error)(-1,
        "BvA.cols=%d does not match A.rows=%d", BvA.M->cols(), A.rows());

    // Apply initial regridding.
    EigenDenseMatrixT B0(*BvA.M * A);        // B0{in}

    // Only apply conservation correction if all of:
    //   a) Matrix is smoothed, so it needs a conservation correction
    //   b) User requested conservation be maintained
    if (BvA.conservative || !force_conservation) {
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

blitz::Array<double,2> WeightedSparse::apply(
    // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
    blitz::Array<double,2> const &A_b,       // A_b{nj} One row per variable
    double fill,    // Fill value for cells not in BvA matrix
    bool force_conservation,
    ibmisc::TmpAlloc &tmp) const
{
    return spsparse::to_blitz<double>(apply_e(A_b, fill), tmp);
}


/** Apply to a single variable */
blitz::Array<double,1> WeightedSparse::apply(
    // WeightedSparse const &BvA,            // BvA_s{ij} smoothed regrid matrix
    blitz::Array<double,1> const &A_b,       // A_b{j} One variable
    double fill,    // Fill value for cells not in BvA matrix
    bool force_conservation,
    ibmisc::TmpAlloc &tmp) const
{
    auto A_b2(ibmisc::reshape<double,1,2>(A_b, {1, A_b.shape()[0]}));
    auto ret2(spsparse::to_blitz(apply_e(A_b2, fill), tmp));
    return ibmisc::reshape<double,2,1>(ret2, {ret2.shape()[1]});
}
// -----------------------------------------------------------------------
void WeightedSparse::ncio(ibmisc::NcIO &ncio,
    std::string const &vname,
    std::array<std::string,2> dim_names)
{

    auto ncdims(ibmisc::get_or_add_dims(ncio,
        {dim_names[0] + ".dense_extent", dim_names[1] + ".dense_extent"},
        {dims[0]->dense_extent(), dims[1]->dense_extent()}));

    // ----------- wM
    std::string matrix_name(dim_names[0] + "v" + dim_names[1]);
    ncio_blitz_alloc<double,1>(ncio, wM, vname+".wM", get_nc_type<double>(),
        {ncdims[0]});

    // --------- M
    ncio_eigen(ncio, *M,
        vname + ".M");

    netCDF::NcVar ncvar = ncio.nc->getVar(vname + ".M.info");
    get_or_put_att(ncvar, ncio.rw,
        "conservative", get_nc_type<bool>(), &conservative, 1);

    // ---- Mw
    ncio_blitz_alloc<double,1>(ncio, Mw, vname + ".Mw", get_nc_type<double>(),
        {ncdims[1]});
}


// -----------------------------------------------------------------
}    // namespace
#endif    // guard
