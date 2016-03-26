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

#include <string>
#include <boost/function.hpp>

class NcFile;


namespace galahad {

    class zd11_f;       // Fortran type, opaque
}

/** @see zdqq_f.f90 */
extern "C" int zd11_put_type_c(galahad::zd11_f &, char const *, int);

namespace galahad {

/**
C++ peer to GALAHAD's derived type zd11_type, which is used to
store sparse matrices.  See galahd namespace documentation for details
on how the peering is achieved.  See GALAHAD QPT documentation for
specifics of how to set the fields in this class.
@see galahad, zd11_type, zd11_x::zd11_c */
class zd11_c {
public :
    zd11_f &main;       ///< Actual storage for this

    int &m;
    int &n;
    int &ne;                ///< Number of non-zero elements

    int * const row;        ///< int[ne]
    int * const col;        ///< int[ne]
    double * const val;     // double[ne]

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
