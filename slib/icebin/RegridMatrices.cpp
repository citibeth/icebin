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
static std::unique_ptr<lintransform::Weighted> compute_AEvI(
    IceRegridder const *regridder,
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &params,
    blitz::Array<double,1> const *elevmaskI,
    UrAE const &AE)
{
printf("BEGIN compute_AEvI scale=%d correctA=%d\n", params.scale, params.correctA);
    std::unique_ptr<lintransform::Eigen> ret(new lintransform::Eigen(dims, true));
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
        std::bind(&IceRegridder::GvI, regridder, _1, elevmaskI),
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

printf("END compute_AEvI\n");
    return ret;
}
// ---------------------------------------------------------
// ---------------------------------------------------------
std::unique_ptr<lintransform::Weighted> compute_IvAE(
    IceRegridder const *regridder,
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &params,
    blitz::Array<double,1> const *elevmaskI,
    UrAE const &AE)
{
printf("BEGIN compute_IvAE\n");
    std::unique_ptr<lintransform::Weighted> ret(new lintransform::Eigen(dims, !params.smooth()));
    SparseSetT * const dimA(ret->dims[1]);    SparseSetT * const dimI(ret->dims[0]);
    SparseSetT _dimG;
    SparseSetT * const dimG(&_dimG);

    if (dimA) dimA->set_sparse_extent(AE.nfull);
    if (dimI) dimI->set_sparse_extent(regridder->nI());
    dimG->set_sparse_extent(regridder->nG());

    // ----- Get the Ur matrices (which determines our dense dimensions)
printf("AA1\n");
    MakeDenseEigenT GvAp_m(
        AE.GvAp,
        {SparsifyTransform::ADD_DENSE},
        {dimG, dimA}, '.');
printf("AA2\n");
    MakeDenseEigenT IvG_m(
        std::bind(&IceRegridder::GvI, regridder, _1, elevmaskI),
        {SparsifyTransform::ADD_DENSE},
        {dimG, dimI}, 'T');
printf("AA3\n");

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
            *dimI, *elevmaskI, ret->wM, params.sigma);
        auto smoothI(to_eigen_sparsematrix(smoothI_t));

        // Smooth the underlying unsmoothed regridding transformation
        ret->M.reset(new EigenSparseMatrixT(smoothI * *ret->M));
    }
printf("END compute_IvAE\n");

    return ret;
}

static std::unique_ptr<lintransform::Weighted> compute_EvA(IceRegridder const *regridder,
    std::array<SparseSetT *,2> dims,
    RegridMatrices::Params const &params, UrAE const &E, UrAE const &A)
{
    std::unique_ptr<lintransform::Weighted> ret(new lintransform::Eigen(dims, true));
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
    blitz::Array<double,1> const &_elevmaskI) const
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
    auto &elevmaskI(rm.tmp.take(blitz::Array<double,1>(_elevmaskI)));

    UrAE urA("UrA", this->nA(),
        std::bind(&IceRegridder::GvAp, regridder, _1, &elevmaskI),
        std::bind(&IceRegridder::sApvA, regridder, _1));

    UrAE urE("UrE", this->nE(),
        std::bind(&IceRegridder::GvEp, regridder, _1, &elevmaskI),
        std::bind(&IceRegridder::sEpvE, regridder, _1));

    // ------- AvI, IvA
    rm.add_regrid("AvI",
        std::bind(&compute_AEvI, regridder, _1, _2, &elevmaskI, urA));
    rm.add_regrid("IvA",
        std::bind(&compute_IvAE, regridder, _1, _2, &elevmaskI, urA));

    // ------- EvI, IvE
    rm.add_regrid("EvI",
        std::bind(&compute_AEvI, regridder, _1, _2, &elevmaskI, urE));
    rm.add_regrid("IvE",
        std::bind(&compute_IvAE, regridder, _1, _2, &elevmaskI, urE));

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
std::unique_ptr<lintransform::Weighted> RegridMatrices::matrix(
    std::string const &spec_name,
    std::array<SparseSetT *,2> dims,
    Params const &params) const
{
    auto &matrix_fn(regrids.at(spec_name));
    auto BvA(matrix_fn(dims, params));
    return BvA;
}
// ----------------------------------------------------------------
// -----------------------------------------------------------------------



// -----------------------------------------------------------------
}    // namespace
#endif    // guard
