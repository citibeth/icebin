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

extern void read_fgrid(
    std::unique_ptr<Grid> &fgridA,    // OUTPUT
    std::string const &fgridA_fname,
    std::string const &fgridA_vname);

extern std::shared_ptr<GCMRegridder_Standard> new_GCMRegridder_Standard(
    Grid const &fgridA,
    std::vector<double> &hcdefs,
    bool _correctA);


/** Instantiates a new C++ object of type GCMRegridder_ModelE.
NOTE: This function is only enabled if BUILD_MODELE is enabled in CMake. */
extern std::shared_ptr<GCMRegridder> new_GCMRegridder_WrapE(
    std::string const &global_ecO,
    std::shared_ptr<GCMRegridder> const &gcmO);

/** Sets the ocean to use for mismatched regridding (icbin::modele::GCMRegridder_ModelE)
@param _gcmA Must be of type GCMRegridder_ModelE *.  Ocean is set for this regridder.
@param foceanAOp_py The ocean (on Ocean grid) as seen by the dynamic ice model.
@param foceanAOm_py The ocean (on Ocean grid) as seen by ModelE. */
extern void GCMRegridder_WrapE_set_focean(
    GCMRegridder *_gcmW,
    PyObject *foceanAOp_py,
    PyObject *foceanAOm_py);

extern PyObject *GCMRegridder_wA(
    GCMRegridder *gcm_regridder,
    std::string const &sheet_name,
    bool native,    // native vs. projected grid
    double fill);

extern void GCMRegridder_add_sheet(GCMRegridder *cself,
    Grid const &fgridA,
    std::string const &name,
    std::string const &gridI_fname, std::string const &gridI_vname,
    std::string const &exgrid_fname, std::string const &exgrid_vname,
    std::string const &sinterp_style);


extern ibmisc::linear::Weighted *RegridMatrices_matrix(RegridMatrices *cself,
    std::string const &spec_name);


PyObject *Hntr_regrid(modele::Hntr const *hntr, PyObject *WTA_py, PyObject *A_py, bool mean_polar);


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
    bool conserve);

/** Allows Python users access to GCMCoupler_Modele::update_topo().
Starting from output of Gary's program (on the Ocean grid), this subroutine
produces a ModelE TOPO file (as internal arrays) on the Atmosphere grid. */
void update_topo(
    // ====== INPUT parameters
    GCMRegridder *_gcmA,
    std::string const &topoO_fname,    // Name of Ocean-based TOPO file (aka Gary)
    PyObject *elevmask_sigmas_py,    // {'greenland' : (elevmaskI<1>, maskI<1>, (sigma_x,signa_y,sigma_z)), ...}
    bool initial_timestep,    // true if this is the first (initialization) timestep
    std::string const &segments,    // [('name', base), ...]
    std::string const &primary_segment,    // [('name', base), ...]
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

