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
#include <icebin/IceRegridder.hpp>
#include <icebin/GCMCoupler.hpp>
#include <icebin/Grid.hpp>
#include <icebin/ElevMask.hpp>
#ifdef BUILD_MODELE
#include <icebin/modele/GCMCoupler_ModelE.hpp>
#endif
#include "icebin_cython.hpp"

using namespace ibmisc;
using namespace ibmisc::cython;
using namespace spsparse;
using namespace icebin;
using namespace icebin::modele;

namespace icebin {
namespace cython {

static double const nan = std::numeric_limits<double>::quiet_NaN();

void read_fgrid(
    std::unique_ptr<Grid> &fgridA,    // Return value goes here
    std::string const &fgridA_fname,
    std::string const &fgridA_vname)
{
    // Read fgridA
    NcIO ncio(fgridA_fname, netCDF::NcFile::read);
    fgridA.reset(new Grid);
    fgridA->ncio(ncio, fgridA_vname);
    ncio.close();
}

std::shared_ptr<GCMRegridder_Standard> new_GCMRegridder_Standard(
    Grid const &fgridA,
    std::vector<double> &hcdefs,
    bool _correctA)
{
    std::shared_ptr<GCMRegridder_Standard> cself(new GCMRegridder_Standard());

    // Put it together
    long nhc = hcdefs.size();
    cself->init(
        AbbrGrid(fgridA),
//        std::move(fgridA),
        std::move(hcdefs),
        Indexing({"A", "HC"}, {0,0}, {(long)fgridA.ndata(), nhc}, {1,0}),
        _correctA);

    return cself;
}

std::shared_ptr<GCMRegridder> new_GCMRegridder_WrapE(
    std::string const &global_ecO,
    std::shared_ptr<GCMRegridder> const &gcmO)
{
#ifdef BUILD_MODELE
    std::shared_ptr<GCMRegridder> gcmW(
        new GCMRegridder_WrapE(
            std::unique_ptr<GCMRegridder_ModelE>(
                new GCMRegridder_ModelE(
                    global_ecO,
                    std::dynamic_pointer_cast<GCMRegridder_Standard>(gcmO)))));
    return gcmW;
#else
    return std::shared_ptr<GCMRegridder>();
#endif
}

void GCMRegridder_WrapE_set_focean(
    GCMRegridder *_gcmW,
    PyObject *foceanAOp_py,
    PyObject *foceanAOm_py)
{
#ifdef BUILD_MODELE
    auto gcmW(dynamic_cast<GCMRegridder_WrapE *>(_gcmW));

    // Check types and convert Numpy Arrays
    size_t nO = gcmW->gcmA->gcmO->nA();
    auto _foceanAOp(np_to_blitz<double,1>(foceanAOp_py, "foceanAOp", {nO}));
    auto _foceanAOm(np_to_blitz<double,1>(foceanAOm_py, "foceanAOm", {nO}));

    // Copy values from Python memory to C++ memory
    gcmW->foceanOp = _foceanAOp;
    gcmW->foceanOm = _foceanAOm;
#endif
}

PyObject *GCMRegridder_wA(
    GCMRegridder *gcm_regridder,
    std::string const &sheet_name,
    bool native,    // native vs. projected grid
    double fill)
{
    PyObject *wA_py = ibmisc::cython::new_pyarray<double,1>(
        std::array<int,1>{gcm_regridder->agridA->dim.sparse_extent()});
    auto wA(np_to_blitz<double,1>(wA_py, "wA", {-1}));
    wA = fill;

    gcm_regridder->wA(
        accum::blitz_existing(wA, DuplicatePolicy::REPLACE), sheet_name, native);

    return wA_py;
}


void GCMRegridder_add_sheet(GCMRegridder *cself,
    Grid const &fgridA,
    std::string const &name,
    std::string const &gridI_fname, std::string const &gridI_vname,
    std::string const &exgrid_fname, std::string const &exgrid_vname,
    std::string const &sinterp_style)
{
    NcIO ncio_I(gridI_fname, netCDF::NcFile::read);
    std::unique_ptr<Grid> fgridI(new Grid);
    fgridI->ncio(ncio_I, "grid");
    ncio_I.close();

    NcIO ncio_exgrid(exgrid_fname, netCDF::NcFile::read);
    std::unique_ptr<Grid> fexgrid(new Grid);
    fexgrid->ncio(ncio_exgrid, exgrid_vname);
    ncio_exgrid.close();

    auto interp_style(parse_enum<InterpStyle>(sinterp_style));

    auto sheet(new_ice_regridder(fgridI->parameterization));
    sheet->init(
        name, *cself->agridA, &fgridA,
        AbbrGrid(*fgridI), ExchangeGrid(*fexgrid),
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


extern linear::Weighted *RegridMatrices_matrix(RegridMatrices *cself, std::string const &spec_name)
{
    return cself->matrix(spec_name).release();
}


// ------------------------------------------------------------
PyObject *Hntr_regrid(Hntr const *hntr, PyObject *WTA_py, PyObject *A_py, bool mean_polar)
{
    // Get Fortran-style (base index = 1) arrays out of this.
    auto WTA(np_to_blitz<double,1>(WTA_py, "WTA", {hntr->Agrid.spec.size()}, blitz::fortranArray));
    auto A(np_to_blitz<double,1>(A_py, "WTA", {hntr->Agrid.spec.size()}, blitz::fortranArray));
    PyObject *B_py(ibmisc::cython::new_pyarray<double,1>(make_array(hntr->Bgrid.spec.size())));
    auto B(np_to_blitz<double,1>(B_py, "WTA", {-1}, blitz::fortranArray));

    hntr->regrid(WTA, A, B, mean_polar);

    return B_py;
}

RegridMatrices *new_regrid_matrices(
    GCMRegridder const *gcm,
    std::string const &sheet_name,
    PyObject *elevmaskI_py,
    // --------- Params
    bool scale,
    bool correctA,
    double sigma_x,
    double sigma_y,
    double sigma_z,
    bool conserve)
{


    auto sheet_index = gcm->ice_regridders().index.at(sheet_name);
    IceRegridder *ice_regridder = &*gcm->ice_regridders()[sheet_index];
    auto elevmaskI(np_to_blitz<double,1>(elevmaskI_py, "elevmaskI", {ice_regridder->nI()}));

    return gcm->regrid_matrices(
        sheet_index, elevmaskI,
        RegridParams(scale, correctA, {sigma_x, sigma_y, sigma_z})).release();
}

std::string to_string(PyObject *str, std::string const &vname)
{
    if (!PyUnicode_Check(str)) (*icebin_error)(-1,
        "Variable %s must be str", vname.c_str());

    Py_ssize_t size;
    char *buf = PyUnicode_AsUTF8AndSize(str, &size);
    return std::string(buf, size);
}

/** Returns: (emI_land, emI_ice) */
PyObject *read_elevmask(std::string const &xfname)
{
    // Read in C++ data structures
    blitz::Array<double,1> emI_land, emI_ice;   // read_elevmask() allocates
    icebin::read_elevmask(xfname, emI_land, emI_ice);

    // Convert to tuple of Numpy arrays
    PyObject *ret = PyTuple_New(2);
        PyTuple_SetItem(ret, 0, copy_blitz_to_np<double,1>(emI_land));
        PyTuple_SetItem(ret, 1, copy_blitz_to_np<double,1>(emI_ice));

    return ret;
}

}}
