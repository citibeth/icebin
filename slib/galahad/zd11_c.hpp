#pragma once

#include <string>
#include <boost/function.hpp>

namespace galahad {

class zd11_f;		// Fortran type, opaque
class zd11_c;		// C type, not opaque to us

}

class NcFile;

extern "C" int zd11_put_type_c(galahad::zd11_f &, char const *, int);

namespace galahad {

// =============== C++ Peer Classes
/** C++ Peer class for the Fortran hsl_zd11_double::zd11_type.
Elements of the Fortran derived type may be safely accessed via the references (scalars) and pointers (arrays) in this class.  See hsl_zd11_double::zd11_type for descriptions of these elements.
<b>NOTE:</b> This class must always match the Fortran derived type zd11_c.
@see hsl_zd11_double::zd11_type, hsl_zd11_double_x::zd11_c */
class zd11_c {
public :
	zd11_f &main;		// Actual storage for this

	int &m;
	int &n;
	int &ne;				// Number of non-zero elements

	/** Array of length ne */
	int * const row;		// int[ne]
	/** Array of length ne */
	int * const col;		// int[ne]
	/** Array of length ne */
	double * const val;		// double[ne]

	/** Construct a dummy instance of this peer class.
	The Fortran subroutine hsl_zd11_double_x::zd11_c_init() is used to fill in
	references to the original Fortran data structure.
	@see zd11_c_init(), hsl_zd11_double::zd11_type */
	zd11_c() : main(*(zd11_f *)0), m(*(int *)0), n(*(int *)0), ne(*(int *)0),
		row(0), col(0), val(0) {}

	/** Set the type parameter in the zd11 data structure.
	@param str Should be 'DENSE', 'COORDINATE', 'SPARSE BY ROWS' or 'DIAGONAL'. */
	int put_type(std::string const &str)
		{ return zd11_put_type_c(main, str.c_str(), str.length()); }

	/** Used to write this data structure to a netCDF file.
	Defines the required variables.  Call the returned boost::function
	later to write the data.
	@param nc NetCDF file to write
	@param vname Variable name to use in writing this sparse matrix.
	@return Function to call later to write the data. */
	boost::function<void()> netcdf_define(NcFile &nc, std::string const &vname);
};

};
