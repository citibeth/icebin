#pragma once

#include <boost/function.hpp>
#include "hsl_zd11_x.hpp"

namespace galahad {

class qtp_problem_f;		// Fortran type, opaque
class qtp_problem_c;			// C++ class

}

class NcFile;

// =============== Fortran Subroutines
extern "C" galahad::qtp_problem_f *qpt_problem_new_c();
extern "C" void qpt_problem_delete_c(galahad::qtp_problem_f *ptr);
extern "C" void qpt_problem_c_init(
		galahad::qtp_problem_c *self, galahad::qtp_problem_f *main,
		int m, int n,
		int A_ne, int H_ne, int eqp_bool);
// =============== C++ Peer Classes

namespace galahad {

/** C++ peer to GALAHAD's derived type galahad_qpt_double::qpt_problem_type. */
class qtp_problem_c {
public :
	qtp_problem_f &main;	/// galahad_qpt_double::qpt_problem_type
	int &m;					/// Number of constraints
	int &n;					/// Number of variables
	double &f;				/// constant term in objective function
	double * const G;		/// double[n] Linear term of objective function
	double * const X_l;		/// double[n] Lower bound on variables
	double * const X_u;		/// double[n] Upper bound on variables
	double * const C_l;		/// double[m] Lower bound, RHS of inequality constraints
	double * const C_u;		/// double[m] Upper bound, RHS of inequality constraints
	double * const C;		/// double[m] RHS of equality constraints
	double * const X;		/// double[n] Value of variables (input & output)
	double * const Y;		/// double[m]
	double * const Z;		/// double[n]

	zd11_c A;				/// m*n Constraints matrix
	zd11_c H;				/// n*n Hessian (quadratic) matrix

	/** Construct new instance.
	Allocates new underlying instance of galahad_qpt_double::qpt_problem_type.  Binds
	to it.  Intializes it.
	@param m Number of constraints.
	@param n Number of variables.
	@param A_ne Number of elements in the constraint matrix.
	@param H_ne Number of elements in the Hessian matrix (the function to be minimized).
	@param eqp Are we preparing for a problem with equality-only constraints? */
	qtp_problem_c(
		int m, int n,
		int A_ne, int H_ne, bool eqp);

	/** Deallocates the underlying istance of galahad_qpt_double::qpt_problem_type. */
	~qtp_problem_c();

	/** Used to write this data structure to a netCDF file.
	Defines the required variables.  Call the returned boost::function
	later to write the data.
	@param nc NetCDF file to write
	@param vname Variable name to use in writing this sparse matrix.
	@return Function to call later to write the data. */
	boost::function<void()> netcdf_define(NcFile &nc, std::string const &vname);
};

}	// namespace galahad
