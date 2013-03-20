#pragma once

#include <memory>
#include <giss/f90blitz.hpp>
#include <giss/SparseMatrix.hpp>
#include <glint2/MatrixMaker.hpp>
#include <glint2/GCMCoupler_MPI.hpp>
#include <glint2/modele/ModelEDomain.hpp>

namespace glint2 {
namespace modele {

struct modele_api {
	std::unique_ptr<MatrixMaker> maker;
	ModelEDomain *domain;	// Points to domain owned by maker

	std::unique_ptr<GCMCoupler_MPI> gcm_coupler;

	// Temporary, until we return the matrix back to the GCM
	std::unique_ptr<giss::VectorSparseMatrix> hp_to_hc;

	// Permanent, used to compute ice sheet SMB
};
}}	// namespace glint2::modele
// ================================================
extern "C" glint2::modele::modele_api *modele_api_new(
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

extern "C" void modele_api_delete(glint2::modele::modele_api *&api);

extern "C"
void modele_api_compute_fhc_c(glint2::modele::modele_api *api,
	giss::F90Array<double, 3> &fhc1h_f,
	giss::F90Array<double, 2> &fgice1_f);

extern "C"
int modele_api_hp_to_hc_part1(glint2::modele::modele_api *api);

extern "C"
void modele_api_hp_to_hc_part2(glint2::modele::modele_api *api,
	giss::F90Array<int, 1> &rows_i_f,
	giss::F90Array<int, 1> &rows_j_f,
	giss::F90Array<int, 1> &rows_k_f,
	giss::F90Array<int, 1> &cols_i_f,
	giss::F90Array<int, 1> &cols_j_f,
	giss::F90Array<int, 1> &cols_k_f,
	giss::F90Array<double, 1> &vals_f);

extern "C"
void modele_api_couple_to_ice(
glint2::modele::modele_api *api,
giss::F90Array<double,3> &smb1hp_f,
giss::F90Array<double,3> &seb1hp_f);
