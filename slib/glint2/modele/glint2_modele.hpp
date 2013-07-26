#pragma once

#include <memory>
#include <giss/f90blitz.hpp>
#include <giss/SparseMatrix.hpp>
#include <glint2/MatrixMaker.hpp>
#include <glint2/GCMCoupler_MPI.hpp>
#include <glint2/modele/ModelEDomain.hpp>

namespace glint2 {
namespace modele {


// ------------------------------------------------------
/** Make a sparse matrix with a vector of theses. */
struct hp_to_ice_rec {
	int row;
	int col_i, col_j, col_k;
	double val;

	hp_to_ice_rec(int _row, int _col_i, int _col_j, int _col_k, double _val) :
		row(_row), col_i(_col_i), col_j(_col_j), col_k(_col_k), val(_val) {}

};
// ------------------------------------------------------

struct glint2_modele {
	std::unique_ptr<MatrixMaker> maker;
	ModelEDomain *domain;	// Points to domain owned by maker

	std::unique_ptr<GCMCoupler_MPI> gcm_coupler;

	std::map<int, std::vector<hp_to_ice_rec>> hp_to_ices;

};
}}	// namespace glint2::modele
// ================================================
extern "C" glint2::modele::glint2_modele *glint2_modele_new(
	char const *maker_fname_f, int maker_fname_len,
	char const *maker_vname_f, int maker_vname_len,

	// Info about the global grid
	int im, int jm,

	// Info about the local grid (C-style indices)
	int i0h, int i1h, int j0h, int j1h,
	int i0, int i1, int j0, int j1,
	int j0s, int j1s,

	// MPI Stuff
	int comm_f, int root,

	// Constants from ModelE
	double LHM, double SHI);

extern "C" void glint2_modele_delete(glint2::modele::glint2_modele *&api);

/** @param replace_fgice_b Should we replace existing fgice1 values with new ones, where the ice sheet overlaps the GCM grid? */
extern "C"
void glint2_modele_compute_fgice_c(glint2::modele::glint2_modele *api,
	int replace_fgice_b,
	giss::F90Array<double, 2> &fgice1_glint2_f,		// OUT
	giss::F90Array<double, 2> &fgice1_f,		// IN/OUT
	giss::F90Array<double, 2> &fgrnd1_f,		// OUT
	giss::F90Array<double, 2> &focean1_f,
	giss::F90Array<double, 2> &flake1_f);

extern "C"
void glint2_modele_init_landice_com_c(glint2::modele::glint2_modele *api,
	giss::F90Array<double, 2> &zatmo1_f,	// IN
	double const BYGRAV,					// IN
	giss::F90Array<double, 2> &fgice1_glint2_f,	// IN
	giss::F90Array<double, 2> &fgice1_f,	// IN
	giss::F90Array<int,3> &used1h_f,		// IN/OUT
	giss::F90Array<double, 3> &fhc1h_f,		// OUT: hp-to-atmosphere
	giss::F90Array<double, 3> &elev1h_f,	// IN/OUT
	int const i0, int const j0, int const i1, int const j1);			// Array bound to write in

extern "C"
void glint2_modele_couple_to_ice_c(
glint2::modele::glint2_modele *api,
int itime,
giss::F90Array<double,3> &smb1hp_f,
giss::F90Array<double,3> &seb1hp_,
giss::F90Array<double,3> &tg21hp_f);
