# this is example-file: configure.ac

# initial information about the project
AC_INIT([glint2],[2.0],[me@myLib.org])
AC_CONFIG_HEADERS([slib/glint2/config.h])
AC_PROG_CXX

m4_include([m4/boost.m4])
#m4_include([m4/ax_python.m4])
#m4_include([m4/ax_lib_cgal_core.m4])
m4_include([m4/ax_prog_doxygen.m4])

# check if the source folder is correct
#AC_CONFIG_SRCDIR([sbin/read_eqp.f90])

# Do shared libraries
#AC_ENABLE_SHARED
#AC_DISABLE_STATIC
#AC_PROG_LIBTOOL(libtool)
LT_INIT([disable-static])


# Checks for programs

# check for C++ compiler and the library compiler
AC_PROG_CXX
#AC_PROG_RANLIB


AC_SUBST([AM_CXXFLAGS], [-frounding-math])

# ----------- Doxygen stuff
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(ON)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(glint2, doxygen.cfg)
# ---------------------------------------------


# automake initialisation (mandatory) and check for minimal automake API version 1.9
AM_INIT_AUTOMAKE([1.9 -Wno-portability])
#AUTOMAKE_OPTIONS = -Wno-portability

# use the C++ compiler for the following checks
AC_LANG([C++])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string])
AC_CHECK_HEADERS([iostream])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# distribute additional compiler and linker flags
# --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
#AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([slib/Makefile])
AC_CONFIG_FILES([sbin/Makefile])
AC_CONFIG_FILES([pyext/Makefile])
AC_CONFIG_FILES([glint2.pc])

# Libraries expected to be available via pkg-config (or otherwise automatically)
#PKG_CHECK_MODULES([NETCDF], [netcdf])
AC_SUBST([NETCDF_CFLAGS])
AC_SUBST([NETCDF_LIBS])

#PKG_CHECK_MODULES([NETCDFF], [netcdf-fortran])
AC_SUBST([NETCDFF_CFLAGS])
AC_SUBST([NETCDFF_LIBS])
AC_SUBST([NETCDFCXX_CFLAGS])
AC_SUBST([NETCDFCXX_LIBS])
#PKG_CHECK_MODULES([PROJ], [proj])
AC_SUBST([PROJ_CFLAGS])
AC_SUBST([PROJ_LIBS])

BOOST_REQUIRE([1.47])
BOOST_THREADS
BOOST_FILESYSTEM
#PKG_CHECK_MODULES([GMP], [gmp])
AC_SUBST([GMP_CFLAGS])
AC_SUBST([GMP_LIBS])
#PKG_CHECK_MODULES([MPFR], [mpfr])
AC_SUBST([MPFR_CFLAGS])
AC_SUBST([MPFR_LIBS])

AC_CHECK_FILE($EIGEN_DIR/Eigen/Core)
AC_SUBST([EIGEN_CFLAGS])

AC_CHECK_FILE($MPI_INCLUDE/mpi.h)
AC_SUBST([MPI_CFLAGS])
AC_SUBST([MPI_LIBS])


#AC_CHECK_FILE($GALAHAD_LIB_DIR/libgalahad.a, \
#	AC_DEFINE(HAVE_GALAHAD,[1],[Indicates presence of GALAHAD library])
#	AM_CONDITIONAL(HAVE_GALAHAD, true)
#	AC_SUBST(GALAHAD_LIB_DIR)
#	AC_SUBST(GALAHAD_INC_DIR), \
#	AC_MSG_ERROR([GALAHAD library not located.])
##	AC_MSG_WARN([GALAHAD library not located.  QP regridding/smoothing will not be enabled.])
#	AM_CONDITIONAL(HAVE_GALAHAD, false))


#AC_CHECK_LIB(mpfr, mpfr_init, , AC_MSG_ERROR([GNU MPFR not found]) )
#AC_CHECK_LIB(gmp, __gmpz_init, , AC_MSG_ERROR([GNU MP not found, see http://gmplib.org]) )
#PKG_CHECK_MODULES([CGAL], [cgal])
AC_SUBST([CGAL_CFLAGS])
AC_SUBST([CGAL_LIBS])

#AX_LIB_CGAL_CORE(, AC_MSG_ERROR([Unable to find CGAL Library]))


# Libraries expected to require manual intervention
#PKG_CHECK_MODULES([BLITZ], [blitz])
AC_SUBST([BLITZ_CFLAGS])
AC_SUBST([BLITZ_LIBS])

AM_PATH_PYTHON
#python
#galahad


AC_PROG_FC


# generate the final Makefile etc.
AC_OUTPUT
