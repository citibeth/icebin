#pragma once

#include <memory>
#include <giss/f90blitz.hpp>
#include <giss/SparseMatrix.hpp>
#include <glint2/MatrixMaker.hpp>
#include <glint2/GCMCoupler_MPI.hpp>
#include <glint2/modele/ModelEDomain.hpp>

namespace glint2 {
namespace modele {


struct glint2_modele_matrix {
	blitz::Array<int, 1> rows_i, rows_j, rows_k;
	blitz::Array<int, 1> cols_i, cols_j, cols_k;
	blitz::Array<double, 1> vals;

	glint2_modele_matrix() {}

	explicit glint2_modele_matrix(int n) :
		rows_i(blitz::Range(1,n)), rows_j(blitz::Range(1,n)), rows_k(blitz::Range(1,n)),
		cols_i(blitz::Range(1,n)), cols_j(blitz::Range(1,n)), cols_k(blitz::Range(1,n)),
		vals(blitz::Range(1,n)) {}


//	explicit glint2_modele_matrix(glint2_modele_matrix_f const &f);
};


struct glint2_modele_matrix_f {
	giss::F90Array<int, 1> rows_i_f;
	giss::F90Array<int, 1> rows_j_f;
	giss::F90Array<int, 1> rows_k_f;
	giss::F90Array<int, 1> cols_i_f;
	giss::F90Array<int, 1> cols_j_f;
	giss::F90Array<int, 1> cols_k_f;
	giss::F90Array<double, 1> vals_f;

	explicit glint2_modele_matrix_f(glint2_modele_matrix &mat) :
		rows_i_f(mat.rows_i), rows_j_f(mat.rows_j), rows_k_f(mat.rows_k),
		cols_i_f(mat.cols_i), cols_j_f(mat.cols_j), cols_k_f(mat.cols_k),
		vals_f(mat.vals) {}

	glint2_modele_matrix to_blitz() {
		glint2_modele_matrix mat;
		mat.rows_i = rows_i_f.to_blitz();
		mat.rows_j = rows_j_f.to_blitz();
		mat.rows_k = rows_k_f.to_blitz();
		mat.cols_i = cols_i_f.to_blitz();
		mat.cols_j = cols_j_f.to_blitz();
		mat.cols_k = cols_k_f.to_blitz();
		mat.vals = vals_f.to_blitz();
		return mat;
	}

};

struct glint2_modele {
	std::unique_ptr<MatrixMaker> maker;
	ModelEDomain *domain;	// Points to domain owned by maker

	std::unique_ptr<GCMCoupler_MPI> gcm_coupler;

	// Temporary, until we return the matrix back to the GCM
	std::unique_ptr<giss::VectorSparseMatrix> hp_to_hc;

	// Permanent, used to compute ice sheet SMB
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
	int comm_f, int root);

extern "C" void glint2_modele_delete(glint2::modele::glint2_modele *&api);

extern "C"
void glint2_modele_compute_fgice_c(glint2::modele::glint2_modele *api,
	giss::F90Array<double, 2> &fgice1_f,		// IN/OUT
	giss::F90Array<double, 2> &fgrnd1_f,		// OUT
	giss::F90Array<double, 2> &focean1_f,
	giss::F90Array<double, 2> &flake1_f);

extern "C"
int glint2_modele_init_landice_com_part1(glint2::modele::glint2_modele *api);

extern "C"
void glint2_modele_init_landice_com_part2(glint2::modele::glint2_modele *api,
	giss::F90Array<double, 3> &fhc1h_f,				// IN/OUT
	giss::F90Array<double, 3> &elevhc_f,			// IN/OUT
	glint2::modele::glint2_modele_matrix_f &hp_to_hc_f,				// OUT
	giss::F90Array<double, 3> &fhp_approx1h_f);		// OUT

extern "C"
void glint2_modele_couple_to_ice(
glint2::modele::glint2_modele *api,
giss::F90Array<double,3> &smb1hp_f,
giss::F90Array<double,3> &seb1hp_f);
