# - Try to find Cgal
# Input Variables
#    CGAL_ROOT
# Once done this will define
#  CGAL_FOUND - System has Cgal
#  CGAL_INCLUDE_DIRS - The Cgal include directories
#  CGAL_LIBRARIES - The libraries needed to use Cgal
##  CGAL_DEFINITIONS - Compiler switches required for using Cgal

find_path(CGAL_ROOT  modules/${CGAL_ARCH}/double/galahad_qpt_double.mod
          HINTS ${CGAL_ROOT})

message(CGAL_ROOT ${CGAL_ROOT})

find_path(CGAL_INCLUDE_DIR CGAL/basic.h
	HINTS ${CGAL_ROOT})

set(CGAL_INCLUDE_DIRS ${CGAL_INCLUDE_DIR})

find_library(CGAL_LIBRARY CGAL
	HINTS ${CGAL_ROOT}/lib)

set(CGAL_LIBRARIES ${CGAL_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CGAL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Cgal  DEFAULT_MSG
                                  CGAL_LIBRARIES CGAL_INCLUDE_DIRS)

mark_as_advanced(CGAL_INCLUDE_DIRS CGAL_LIBRARIES )




#
## Find the CGAL includes and client library
## This module defines
##  CGAL_INCLUDE_DIR, where to find CGAL.h
##  CGAL_LIBRARIES, the libraries needed to use CGAL.
##  CGAL_FOUND, If false, do not try to use CGAL.
#
#if(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
#   set(CGAL_FOUND TRUE)
#
#else(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
#
# FIND_PATH(CGAL_INCLUDE_DIR CGAL/basic.h
#      /usr/include
#      /usr/local/include
#      $ENV{ProgramFiles}/CGAL/*/include
#      $ENV{SystemDrive}/CGAL/*/include
#      )
#
#  find_library(CGAL_LIBRARIES NAMES CGAL libCGAL
#     PATHS
#     /usr/lib
#     /usr/local/lib
#     /usr/lib/CGAL
#     /usr/lib64
#     /usr/local/lib64
#     /usr/lib64/CGAL
#     $ENV{ProgramFiles}/CGAL/*/lib
#     $ENV{SystemDrive}/CGAL/*/lib
#     )
#     
#  if(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
#    set(CGAL_FOUND TRUE)
#    message(STATUS "Found CGAL: ${CGAL_INCLUDE_DIR}, ${CGAL_LIBRARIES}")
#    INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIR} $ENV{CGAL_CFG})
#  else(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
#    set(CGAL_FOUND FALSE)
#    message(STATUS "CGAL not found.")
#  endif(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
#
#  mark_as_advanced(CGAL_INCLUDE_DIR CGAL_LIBRARIES)
#
#endif(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
#