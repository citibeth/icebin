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
#include "icebin_cython.hpp"

using namespace ibmisc;
using namespace ibmisc::cython;
using namespace spsparse;

namespace icebin {
namespace cython {

void GCMRegridder_init(GCMRegridder *cself,
    std::string const &gridA_fname,
    std::string const &gridA_vname,
    std::vector<double> &hcdefs,
    bool _correctA)
{
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

}


void GCMRegridder_add_sheet(GCMRegridder *cself,
    std::string name,
    std::string const &gridI_fname, std::string const &gridI_vname,
    std::string const &exgrid_fname, std::string const &exgrid_vname,
    std::string const &sinterp_style,
    PyObject *elevI_py)
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
    auto elevI(np_to_blitz<double,1>(elevI_py, "elevI", {gridI->ndata()}));

    auto sheet(new_ice_regridder(gridI->parameterization));
    sheet->init(name, std::move(gridI), std::move(exgrid),
        interp_style, elevI);
    cself->add_sheet(std::move(sheet));
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

PyObject *RegridMatrices_regrid(RegridMatrices *cself, std::string const &spec_name, bool scale, bool correctA, double sigma)
{
    std::array<SparseSetT,2> dims;
    std::unique_ptr<WeightedSparse> Mw(cself->regrid(spec_name,
        {&dims[0], &dims[1]},
        RegridMatrices::Params(scale, correctA, sigma)));

    // ----- Convert a dense vector w/ dense indices to a dense vector with sparse indices
    // Allocate the output Python vector
    PyObject *weight_py = ibmisc::cython::new_pyarray<double,1>(
        std::array<npy_intp,1>{dims[0].sparse_extent()});
    // Get a Blitz view on it
    auto weight_pyb(np_to_blitz<double,1>(weight_py, "weight_py", {-1}));

    // Copy, while translating the dimension
    spcopy(
        accum::to_sparse(make_array(&dims[0]),
        accum::blitz_existing(weight_pyb)),
        Mw->weight);

    // Convert a sparse matrix w/ dense indices to a sparse matrix with sparse indices
    TupleListT<2> Mw_sp;
    spcopy(
        accum::to_sparse(make_array(&dims[0], &dims[1]),
        accum::ref(Mw_sp)),
        *Mw->M);
    PyObject *M_py = ibmisc::cython::spsparse_to_tuple(Mw_sp);
    fflush(stdout);
    return Py_BuildValue("OO", M_py, weight_py);
}


}}

