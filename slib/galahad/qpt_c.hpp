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

#include <boost/function.hpp>
#include <galahad/zd11_c.hpp>

class NcFile;

namespace galahad {
	class qpt_problem_f;		// Fortran type, opaque
	class qpt_problem_c;		// Forward declaration for functions below
}

/** @see qpt_f.f90 */
extern "C" void qpt_problem_alloc_h(
	galahad::qpt_problem_c *this_c, int H_ne);

/** @see qpt_f.f90 */
extern "C" void qpt_problem_alloc_a(
	galahad::qpt_problem_c *this_c, int A_ne);


namespace galahad {

/** C++ peer to GALAHAD's derived type qpt_problem_type (in GALAHAD documentation).  See
galahd namespace documentation for details on how the peering is
achieved.
@see galahad, qpt_problem_type, qpt_x::qpt_problem_c */

class qpt_problem_c {
public :
	qpt_problem_f &this_f;	//!< galahad_qpt_double::qpt_problem_type
	int &m;					//!< Number of constraints
	int &n;					//!< Number of variables
	double &f;				//!< constant term in objective function
	int &Hessian_kind;		//!< See QPT documenation

	/** Description of G */
	double * const G;		// double[n] Linear term of objective function
	double * const X_l;		//!< double[n] Lower bound on variables
	double * const X_u;		//!< double[n] Upper bound on variables
	double * const C_l;		//!< double[m] Lower bound, RHS of inequality constraints
	double * const C_u;		//!< double[m] Upper bound, RHS of inequality constraints
	double * const C;		//!< double[m] RHS of equality constraints
	double * const X;		//!< double[n] Value of variables (input & output)
	double * const Y;		//!< double[m]
	double * const Z;		//!< double[n]
//	double * const WEIGHT;
	double * const X0;		//!< double[n] For weighted least square distance problem

	zd11_c A;				//!< m*n Constraints matrix
	zd11_c H;				//!< n*n Hessian (quadratic) matrix

	/** Construct new instance.
	Allocates new underlying Fortran instance of QPT_problem_type.
	    Binds to it.  Intializes it.  Note that this does \emph{not}
	    allocate the H and A matrices; that must be done separately
	    through alloc_H() and alloc_A().
	@param m Number of constraints.
	@param n Number of variables.
	@param A_ne Number of elements in the constraint matrix.
	@param H_ne Number of elements in the Hessian matrix (the function to be minimized).
	@param eqp Are we preparing for a problem with equality-only constraints?
	@param Hessian_kind See the GALAHAD documentation for QPT.  NOTE:
        Hessian_kind is not used in the GALAHAD EQP module, do not try
        to change it away from -1 in this case. */
	qpt_problem_c(
		int m, int n,
		bool eqp,
		int Hessian_kind = -1);

	/** Calls Fortran code to allocate space for the matrix H.  This
	allows late binding of the number of elements to be stored in that
	matrix.

	@see qpt_x::qpt_problem_alloc_h */
	void alloc_H(int H_ne)
		{ qpt_problem_alloc_h(this, H_ne); }

	/** Calls Fortran code to allocate space for the matrix A.  This
	allows late binding of the number of elements to be stored in that
	matrix.

	@see qpt_x::qpt_problem_alloc_a */
	void alloc_A(int A_ne)
		{ qpt_problem_alloc_a(this, A_ne); }


	/** Deallocates the underlying istance of qpt_x::qpt_problem_type. */
	~qpt_problem_c();

	/** Evaluates the objective function, apart from GALAHAD.
	@param x x[n] Pointer to array holding vector value at which
	objective function is to be evaluated. */
	double eval_objective(double const *x);

	/** Used to write this data structure to a netCDF file.
	Defines the required variables.  Call the returned boost::function
	later to write the data.
	@param nc NetCDF file to write
	@param vname Variable name to use in writing this sparse matrix.
	@return Function to call later to write the data. */
	boost::function<void()> netcdf_define(NcFile &nc, std::string const &vname);
};

}	// namespace galahad
