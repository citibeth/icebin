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

#pragma once

#include <Python.h>
#include <ibmisc/cython.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/modele/hntr.hpp>

namespace icebin {
namespace cython {

extern std::shared_ptr<GCMRegridder_Standard> new_GCMRegridder_Standard(
    std::string const &gridA_fname,
    std::string const &gridA_vname,
    std::vector<double> &hcdefs,
    bool _correctA);

extern std::shared_ptr<GCMRegridder> new_GCMRegridder_ModelE(
    std::shared_ptr<GCMRegridder> const &gcmO,
    PyObject *foceanAOp_py,
    PyObject *foceanAOm_py);

extern PyObject *GCMRegridder_wA(
    GCMRegridder *gcm_regridder,
    std::string const &sheet_name,
    bool native,    // native vs. projected grid
    double fill);

extern void GCMRegridder_add_sheet(GCMRegridder *cself,
    std::string const &name,
    std::string const &gridI_fname, std::string const &gridI_vname,
    std::string const &exgrid_fname, std::string const &exgrid_vname,
    std::string const &sinterp_style);


/** Wraps WeightedSparse to keep around the dense/sparse dimension
    translators */
struct CythonWeightedSparse {
    std::array<SparseSetT,2> dims;
    std::unique_ptr<WeightedSparse> RM;

    /** @return 2-D shape of the stored matrix */
    PyObject *shape();
};


extern CythonWeightedSparse *RegridMatrices_matrix(RegridMatrices *cself,
    std::string const &spec_name, bool scale, bool correctA,
    double sigma_x, double sigma_y, double sigma_z, bool conserve);

extern PyObject *CythonWeightedSparse_apply(
    CythonWeightedSparse *BvA,
    PyObject *A_s_py,
    double fill, bool force_conservation);            // A_b{nj_s} One row per variable


PyObject *CythonWeightedSparse_dense_extent(CythonWeightedSparse const *cself);
PyObject *CythonWeightedSparse_sparse_extent(CythonWeightedSparse const *cself);

PyObject *CythonWeightedSparse_to_tuple(CythonWeightedSparse *cself);

void coo_matvec(PyObject *yy_py, PyObject *xx_py, bool ignore_nan,
    size_t M_nrow, size_t M_ncol, PyObject *M_row_py, PyObject *M_col_py, PyObject *M_data_py);

PyObject *Hntr_regrid(modele::Hntr const *hntr, PyObject *WTA_py, PyObject *A_py, bool mean_polar);


RegridMatrices *new_regrid_matrices(GCMRegridder const *gcm, std::string const &sheet_name, PyObject *elevI_py);

void update_topo(
    // ====== INPUT parameters
    GCMRegridder *_gcmA,
    std::string const &topoO_fname,    // Name of Ocean-based TOPO file (aka Gary)
    PyObject *elevmask_sigmas_py,    // {'greenland' : (elevI<1>, maskI<1>, (sigma_x,signa_y,sigma_z)), ...}
    bool initial_timestep,    // true if this is the first (initialization) timestep
    std::string const &segments,    // [('name', base), ...]
    // ===== OUTPUT parameters (variables come from GCMCoupler); must be pre-allocated
    PyObject *fhc_py,         // blitz::Array<double,3> fhc;
    PyObject *underice_py,    // blitz::Array<int,3> underice;
    PyObject *elevE_py,       // blitz::Array<double,3> elevE;
    // i,j arrays on Atmosphere grid
    PyObject *focean_py,      // blitz::Array<double,2> focean;
    PyObject *flake_py,       // blitz::Array<double,2> flake;
    PyObject *fgrnd_py,       // blitz::Array<double,2> fgrnd;    // Alt: fearth0
    PyObject *fgice_py,       // blitz::Array<double,2> fgice;    // Alt: flice
    PyObject *zatmo_py,       // blitz::Array<double,2> zatmo;      // i,j

    PyObject *foceanOm0_py);   // blitz::Array<double,1> foceanOm0,

}}

