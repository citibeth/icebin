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

namespace icebin {
namespace cython {

extern void GCMRegridder_init(GCMRegridder *cself,
    std::string const &gridA_fname,
    std::string const &gridA_vname,
    std::vector<double> &hpdefs,
    bool _correctA);

extern void GCMRegridder_add_sheet(GCMRegridder *cself,
    std::string name,
    std::string const &gridI_fname, std::string const &gridI_vname,
    std::string const &exgrid_fname, std::string const &exgrid_vname,
    std::string const &sinterp_style,
    PyObject *elevI_py, PyObject *maskI_py);

extern PyObject *RegridMatrices_regrid(RegridMatrices *cself, std::string const &spec_name, bool scale, bool correctA, double sigma);

void coo_matvec(PyObject *yy_py, PyObject *xx_py, bool ignore_nan,
    size_t M_nrow, size_t M_ncol, PyObject *M_row_py, PyObject *M_col_py, PyObject *M_data_py);


}}

