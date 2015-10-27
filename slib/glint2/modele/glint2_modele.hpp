/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <memory>
#include <giss/f90blitz.hpp>
#include <giss/SparseMatrix.hpp>
#include <glint2/MatrixMaker.hpp>
#include <glint2/GCMCoupler.hpp>
#include <glint2/modele/GCMCoupler_ModelE.hpp>
#include <glint2/modele/ModelEDomain.hpp>

namespace glint2 {
namespace modele {

// The ModelE interface to GLINT2.  This is called directly by ModelE.
// ------------------------------------------------------

/** Make a sparse matrix with a vector of these.  This struct is
specific to ModelE's (i,j,k) indexing scheme. */
struct hp_to_ice_rec {
	int row;
	int col_i, col_j, col_k;
	double val;

	hp_to_ice_rec(int _row, int _col_i, int _col_j, int _col_k, double _val) :
		row(_row), col_i(_col_i), col_j(_col_j), col_k(_col_k), val(_val) {}

};
// ------------------------------------------------------

struct glint2_modele {
	ModelEDomain *domain;	// Points to domain owned by maker

	double dtsrc;			// Size of ModelE timestep
	GCMCoupler_ModelE gcm_coupler;

	/** Vectors handed to us by ModelE, before we pass them along to
	Glint2 when coupling. */
	std::vector<blitz::Array<double, 3>> gcm_outputs;

	/** Position with the gcm_inputs array (passed from GCM)
	where each variable in the gcm_inputs contract starts.
	This has a sentinel on the end, thus indicating the extent of
	each variable as well. */
	std::vector<int> gcm_inputs_ihp;

	/** The matrix used for each IceModel, used to convert from
	the elevation grid to ice grid.  Each hp_to_ice_rec is one
	non-zero element of the matrix. */
	std::map<int, std::vector<hp_to_ice_rec>> hp_to_ices;

	/** Last time the coupler was called (or start of run) */
	int itime_last;

	glint2_modele()
	{
		/** Sentinel */
		gcm_inputs_ihp.push_back(0);
	}
};

}}	// namespace glint2::modele
// ================================================
/** Just allocate */
extern "C" glint2::modele::glint2_modele *new_glint2_modele_c();

/** Set one constant from ModelE into GLINT2.  This is used as a callback
function for ModelE's set_constant() Fortran subroutine. */
extern "C" void glint2_modele_set_const(
	glint2::modele::glint2_modele *api,
	char const *name_f, int name_len,
	double val,
	char const *units_f, int units_len,
	char const *description_f, int description_len);

/** Tell the GCM how many elevation points are involved in this
configuration, INCLUDING the "legacy" elevation point, which is
used by ModelE but not Glint2. */
extern "C"
int glint2_modele_nhp_gcm(glint2::modele::glint2_modele const *api);

/** Inform Glint2 about a Fortran variable used to hold inputs to the
GCM (regridded from the ice model output).  This is called from
ModelE, as a way to ensure that ModelE and Glint2 are working from
the same set of variables.
*/
extern "C"
int glint2_modele_add_gcm_input(
glint2::modele::glint2_modele *api,
char const *field_name_f, int field_name_len,
char const *units_f, int units_len,
char const *grid_f, int grid_len,
int initial,	// bool
char const *long_name_f, int long_name_len);


// First init to be called after new
extern "C" void glint2_modele_init0(
	glint2::modele::glint2_modele *api,
	char const *maker_fname_f, int maker_fname_len,
	char const *maker_vname_f, int maker_vname_len,

	// Info about the global grid
	int im, int jm,

	// Info about the local grid (C-style indices)
	int i0h, int i1h, int j0h, int j1h,
	int i0, int i1, int j0, int j1,
	int j0s, int j1s,

	// MPI Stuff
	MPI_Fint comm_f, int root,

	// API  control
	int write_constants);

extern "C" void glint2_modele_delete(glint2::modele::glint2_modele *&api);

extern "C"
void glint2_modele_set_start_time(glint2::modele::glint2_modele *api,
	int iyear1, int itimei, double dtsrc);

extern "C"
void glint2_modele_get_flice_im_c(glint2::modele::glint2_modele *api,
	giss::F90Array<double, 2> &flice1_im_f);		// OUT

/** Produces the (dense) FHC_IM array from the (sparse) hp_to_atm
coming from raw Glint2. */
extern "C"
void glint2_modele_get_fhc_im_c(glint2::modele::glint2_modele *api,
	giss::F90Array<double, 3> &fhc_im1h_f);	// OUT

extern "C"
void glint2_modele_get_elevhp_im_c(glint2::modele::glint2_modele *api,
	giss::F90Array<double, 3> &elev1h_f);	// IN/OUT

extern "C"
void glint2_modele_init_hp_to_ices(glint2::modele::glint2_modele *api);

extern "C"
void glint2_modele_couple_to_ice_c(
glint2::modele::glint2_modele *api,
int itime,			// ModelE itime counter
giss::F90Array<double,3> &gcm_inputs_d_f);

extern "C"
void  glint2_modele_get_initial_state_c(
glint2::modele::glint2_modele *api,
giss::F90Array<double,3> &gcm_inputs_d_f);
