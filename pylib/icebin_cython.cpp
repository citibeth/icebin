/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include <algorithm>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/cython.hpp>
#include <spsparse/SparseSet.hpp>
#ifdef BUILD_MODELE
#include <icebin/modele/GCMRegridder_ModelE.hpp>
#endif
#include "icebin_cython.hpp"

using namespace ibmisc;
using namespace ibmisc::cython;
using namespace spsparse;
using namespace icebin::modele;

namespace icebin {
namespace cython {

static double const nan = std::numeric_limits<double>::quiet_NaN();

std::shared_ptr<GCMRegridder_Standard> new_GCMRegridder_Standard(
    std::string const &gridA_fname,
    std::string const &gridA_vname,
    std::vector<double> &hcdefs,
    bool _correctA)
{
    std::shared_ptr<GCMRegridder_Standard> cself(new GCMRegridder_Standard());

    // Read gridA
    NcIO ncio(gridA_fname, netCDF::NcFile::read);
    std::unique_ptr<Grid> gridA = new_grid(ncio, gridA_vname);
    gridA->ncio(ncio, gridA_vname);
    ncio.close();

    // Put it together
    long nhc = hcdefs.size();
    cself->init(
        std::move(gridA),
        std::move(hcdefs),
        Indexing({"A", "HC"}, {0,0}, {gridA->ndata(), nhc}, {1,0}),
        _correctA);

    return cself;
}

std::shared_ptr<GCMRegridder> new_GCMRegridder_ModelE(
    std::shared_ptr<GCMRegridder> const &gcmO,
    PyObject *foceanAOp_py,
    PyObject *foceanAOm_py)
{
#ifdef BUILD_MODELE
    std::shared_ptr<GCMRegridder_ModelE> gcmA(new GCMRegridder_ModelE(gcmO));

    // Check types and convert Numpy Arrays
    size_t nO = gcmO->nA();
    auto _foceanAOp(np_to_blitz<double,1>(foceanAOp_py, "foceanAOp", {nO}));
    auto _foceanAOm(np_to_blitz<double,1>(foceanAOm_py, "foceanAOm", {nO}));

    // Copy values from Python memory to C++ memory
    blitz::Array<double,1> foceanAOp(_foceanAOp.shape());
    foceanAOp = _foceanAOp;
    blitz::Array<double,1> foceanAOm(_foceanAOm.shape());
    foceanAOm = _foceanAOm;

    gcmA->set_focean(foceanAOp, foceanAOm);
    return gcmA;
#else
    return std::shared_ptr<GCMRegridder_ModelE>();
#endif
}

PyObject *GCMRegridder_wA(
    GCMRegridder *gcm_regridder,
    std::string const &sheet_name,
    bool native,    // native vs. projected grid
    double fill)
{
    PyObject *wA_py = ibmisc::cython::new_pyarray<double,1>(
        std::array<int,1>{gcm_regridder->gridA->ndata()});
    auto wA(np_to_blitz<double,1>(wA_py, "wA", {-1}));
    wA = fill;

    gcm_regridder->wA(
        accum::blitz_existing(wA, DuplicatePolicy::REPLACE), sheet_name, native);

    return wA_py;
}


void GCMRegridder_add_sheet(GCMRegridder *cself,
    std::string const &name,
    std::string const &gridI_fname, std::string const &gridI_vname,
    std::string const &exgrid_fname, std::string const &exgrid_vname,
    std::string const &sinterp_style)
{
    NcIO ncio_I(gridI_fname, netCDF::NcFile::read);
    std::unique_ptr<Grid> gridI(new_grid(ncio_I, "grid"));
    gridI->ncio(ncio_I, "grid");
    ncio_I.close();

    NcIO ncio_exgrid(exgrid_fname, netCDF::NcFile::read);
    std::unique_ptr<Grid> exgrid(new_grid(ncio_exgrid, exgrid_vname));
    exgrid->ncio(ncio_exgrid, exgrid_vname);
    ncio_exgrid.close();

    auto interp_style(parse_enum<InterpStyle>(sinterp_style));

    auto sheet(new_ice_regridder(gridI->parameterization));
    sheet->init(name, std::move(gridI), std::move(exgrid),
        interp_style);
    dynamic_cast<GCMRegridder_Standard *>(cself)
        ->add_sheet(std::move(sheet));
}

/** Computes yy = M xx.  yy is allocated, not necessarily set. */
void coo_matvec(PyObject *yy_py, PyObject *xx_py, bool ignore_nan,
    size_t M_nrow, size_t M_ncol, PyObject *M_row_py, PyObject *M_col_py, PyObject *M_data_py)
{
    auto xx(np_to_blitz<double,1>(xx_py, "xx", {M_ncol}));
    auto yy(np_to_blitz<double,1>(yy_py, "yy", {M_nrow}));
    auto M_row(np_to_blitz<int,1>(M_row_py, "M_row_py", {-1}));
    auto M_col(np_to_blitz<int,1>(M_col_py, "M_col_py", {-1}));
    auto M_data(np_to_blitz<double,1>(M_data_py, "M_data_py", {-1}));

    // Keep track of which items we've written to.
    std::vector<bool> written(M_nrow, false);

    // Do the multiplication, and we're done!
    int nnz = M_data.size();
    for (int n=0; n<nnz; ++n) {
        size_t row = M_row(n);
        size_t col = M_col(n);
        double data = M_data(n);

        // Ignore NaN in input vector
        if (ignore_nan && std::isnan(xx(col))) continue;

        // Just do Snowdrift-style "REPLACE".  "MERGE" was never used.
        double old_yy;
        if (written[row]) {
            old_yy = yy(row);
        } else {
            old_yy = 0;
            written[row] = true;
        }
        yy(row) = old_yy + data * xx(col);
    }

}


PyObject *CythonWeightedSparse::shape()
{
    return py_build_tuple<pytype_long, long, 2>(
        std::array<long,2>{
            RM->dims[0]->sparse_extent(),
            RM->dims[1]->sparse_extent()
        }
    );
}


extern CythonWeightedSparse *RegridMatrices_matrix(RegridMatrices *cself, std::string const &spec_name, bool scale, bool correctA, double sigma_x, double sigma_y, double sigma_z, bool conserve)
{
    std::unique_ptr<CythonWeightedSparse> CRM(new CythonWeightedSparse());

    CRM->RM = cself->matrix(spec_name,
        {&CRM->dims[0], &CRM->dims[1]},
        RegridMatrices::Params(scale, correctA, {sigma_x, sigma_y, sigma_z}));

    CythonWeightedSparse *ret = CRM.release();
    return ret;
}

extern PyObject *CythonWeightedSparse_apply(
    CythonWeightedSparse *BvA,
    PyObject *A_s_py,            // A_b{nj_s} One row per variable
    double fill, bool force_conservation)
{
    // |j_s| = size of sparse input vector space (A_s)
    // |j_d] = size of dense input vector space (A_d)
    // |n| = number of variables being processed together

    // Allocate dense A matrix
    auto &bdim(BvA->dims[0]);
    auto &adim(BvA->dims[1]);
    auto A_s(np_to_blitz<double,2>(A_s_py, "A", {-1,-1}));    // Sparse indexed dense vector
    int n_n = A_s.extent(0);

    // Densify the A matrix
    blitz::Array<double,2> A_d(n_n, adim.dense_extent());
    for (int j_d=0; j_d < adim.dense_extent(); ++j_d) {
        int j_s = adim.to_sparse(j_d);
        for (int n=0; n < n_n; ++n) {
            A_d(n,j_d) = A_s(n,j_s);
        }
    }

    // Apply...
    auto B_d_eigen(BvA->RM->apply_e(A_d, fill, force_conservation));    // Column major indexing

    // Allocate output vector and get a Blitz view
    // We will copy from the Eigen data structure to the Python
    PyObject *B_s_py = ibmisc::cython::new_pyarray<double,2>(
        std::array<int,2>{n_n, bdim.sparse_extent()});
    auto B_s(np_to_blitz<double,2>(B_s_py, "B_s_py", {-1,-1}));

    // Sparsify the output B
    for (int n=0; n < n_n; ++n)
    for (int j_s=0; j_s < bdim.sparse_extent(); ++j_s) {
        B_s(n,j_s) = fill;
    }

    for (int j_d=0; j_d < bdim.dense_extent(); ++j_d) {
        int j_s = bdim.to_sparse(j_d);
        for (int n=0; n < n_n; ++n) {
            B_s(n,j_s) = B_d_eigen(j_d,n);
        }
    }

    return B_s_py;
}

PyObject *CythonWeightedSparse_dense_extent(CythonWeightedSparse const *cself)
{
    return Py_BuildValue("ll",
        (long)cself->dims[0].dense_extent(),
        (long)cself->dims[1].dense_extent());
}

PyObject *CythonWeightedSparse_sparse_extent(CythonWeightedSparse const *cself)
{
    return Py_BuildValue("ll",
        (long)cself->dims[0].sparse_extent(),
        (long)cself->dims[1].sparse_extent());
}


PyObject *CythonWeightedSparse_to_tuple(CythonWeightedSparse *cself)
{
    // ----- Convert a dense vector w/ dense indices to a dense vector with sparse indices
    if (cself->dims[0].sparse_extent() <= 0) (*icebin_error)(-1,
        "dims[0].sparse_extent() must be > 0");

    // ---------- wM
    // Allocate the output Python vector
    PyObject *wM_py = ibmisc::cython::new_pyarray<double,1>(
        std::array<int,1>{cself->dims[0].sparse_extent()});
    // Get a Blitz view on it
    auto wM_pyb(np_to_blitz<double,1>(wM_py, "wM_py", {-1}));
    // Copy, while translating the dimension
    spcopy(
        accum::to_sparse(make_array(&cself->dims[0]),
        accum::blitz_existing(wM_pyb)),
        cself->RM->wM);

    // -------------- M
    // Convert a sparse matrix w/ dense indices to a sparse matrix with sparse indices
    TupleListT<2> RM_sp;
    spcopy(
        accum::to_sparse(make_array(&cself->dims[0], &cself->dims[1]),
        accum::ref(RM_sp)),
        *cself->RM->M);
    PyObject *M_py = ibmisc::cython::spsparse_to_tuple(RM_sp);

    // ---------- Mw
    // Allocate the output Python vector
    PyObject *Mw_py = ibmisc::cython::new_pyarray<double,1>(
        std::array<int,1>{cself->dims[1].sparse_extent()});
    // Get a Blitz view on it
    auto Mw_pyb(np_to_blitz<double,1>(Mw_py, "Mw_py", {-1}));
    // Copy, while translating the dimension
    spcopy(
        accum::to_sparse(make_array(&cself->dims[1]),
        accum::blitz_existing(Mw_pyb)),
        cself->RM->Mw);

    return Py_BuildValue("OOO", wM_py, M_py, Mw_py);
}

// ------------------------------------------------------------
PyObject *Hntr_regrid(Hntr const *hntr, PyObject *WTA_py, PyObject *A_py, bool mean_polar)
{
    // Get Fortran-style (base index = 1) arrays out of this.
    auto WTA(np_to_blitz<double,1>(WTA_py, "WTA", {hntr->Agrid.size()}, blitz::fortranArray));
    auto A(np_to_blitz<double,1>(A_py, "WTA", {hntr->Agrid.size()}, blitz::fortranArray));
    PyObject *B_py(ibmisc::cython::new_pyarray<double,1>(make_array(hntr->Bgrid.size())));
    auto B(np_to_blitz<double,1>(B_py, "WTA", {-1}, blitz::fortranArray));

    hntr->regrid(WTA, A, B, mean_polar);

    return B_py;
}

RegridMatrices *new_regrid_matrices(GCMRegridder const *gcm, std::string const &sheet_name, PyObject *elevI_py)
{
    IceRegridder *ice_regridder = gcm->ice_regridder(sheet_name);
    auto elevI(np_to_blitz<double,1>(elevI_py, "elevI", {ice_regridder->nI()}));
    return new RegridMatrices(gcm->regrid_matrices(sheet_name, &elevI));
}


}}
