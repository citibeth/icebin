# This file contains CMake macros used in the root CMakeLists.txt

macro(glint2_find_prerequisites)
  # PETSc
  find_package (PETSc)
  if (NOT PETSC_FOUND)
    get_filename_component(pcc ${PETSC_COMPILER} REALPATH)
    get_filename_component(cc ${CMAKE_C_COMPILER} REALPATH)
    if (NOT ${pcc} STREQUAL ${cc})
      message(WARNING
        "PETSC_COMPILER does not match CMAKE_C_COMPILER\n"
	"  PETSC_COMPILER=${PETSC_COMPILER}\n"
	"  CMAKE_C_COMPILER=${CMAKE_C_COMPILER}\n"
	"Try running \n"
	"  rm CMakeCache.txt && cmake -DCMAKE_C_COMPILER=${PETSC_COMPILER} ${CMAKE_SOURCE_DIR}")
    endif()
    message(FATAL_ERROR  "GLINT2 configuration failed: PETSc was not found.")
  endif()

  if ((DEFINED PETSC_VERSION) AND (PETSC_VERSION VERSION_LESS 3.3))
    # Force GLINT2 to look for PETSc again if the version we just found
    # is too old:
    set(PETSC_CURRENT "OFF" CACHE BOOL "" FORCE)
    # Stop with an error message.
    message(FATAL_ERROR "\nGLINT2 requires PETSc version 3.3 or newer (found ${PETSC_VERSION}).\n\n")
  endif()

  # MPI
  # Use the PETSc compiler as a hint when looking for an MPI compiler
  # FindMPI.cmake changed between 2.8.4 and 2.8.5, so we try to support both...
  if (${CMAKE_VERSION} VERSION_LESS "2.8.5")
    set (MPI_COMPILER ${PETSC_COMPILER} CACHE FILEPATH "MPI compiler. Used only to detect MPI compilation flags.")
    find_package (MPI REQUIRED)

    set (MPI_C_INCLUDE_PATH "${MPI_INCLUDE_PATH}" CACHE STRING "MPI include directories (semicolon-separated list)")
    set (MPI_C_LIBRARIES "${MPI_LIBRARY};${MPI_EXTRA_LIBRARY}" CACHE STRING "MPI libraries (semicolon-separated list)")
    mark_as_advanced(MPI_C_INCLUDE_PATH MPI_C_LIBRARIES)
    message (STATUS
      "Note: Please upgrade CMake to version 2.8.5 or later if the build fails with undefined symbols related to MPI.")
  else ()
    set (MPI_C_COMPILER ${PETSC_COMPILER} CACHE FILEPATH "MPI compiler. Used only to detect MPI compilation flags.")
    find_package (MPI REQUIRED)
  endif()

  # Other required libraries
  find_package (UDUNITS2 REQUIRED)
  find_package (GSL REQUIRED)
  find_package (NetCDF REQUIRED)
  find_package (CGAL REQUIRED)

message("CGAL " ${CGAL_INCLUDE_DIR})
message("CGAL " ${CGAL_LIBRARIES})


  # Optional libraries
#  find_package (PNetCDF)
  find_package (FFTW REQUIRED)
  find_package (PROJ4)
#  find_package (TAO)
#  # Try to find netcdf_par.h. We assume that NetCDF was compiled with
#  # parallel I/O if this header is present.
#  find_file(NETCDF_PAR_H netcdf_par.h HINTS ${NETCDF_INCLUDES} NO_DEFAULT_PATH)

#  # Set default values for build options
#  if (NOT NETCDF_PAR_H)
#    set (Glint2_USE_PARALLEL_NETCDF4 OFF CACHE BOOL "Enables parallel NetCDF-4 I/O." FORCE)
#    message(STATUS "Selected NetCDF library does not support parallel I/O.")
#  endif()
#
#  if (NOT PNETCDF_FOUND)
#    set (Glint2_USE_PNETCDF OFF CACHE BOOL "Enables parallel NetCDF-3 I/O using PnetCDF." FORCE)
#  endif()
#
#  if (NOT PROJ4_FOUND)
#    set (Glint2_USE_PROJ4 OFF CACHE BOOL "Use Proj.4 to compute cell areas, longitude, and latitude." FORCE)
#  endif()
#
#  if (NOT TAO_FOUND)
#    set (Glint2_USE_TAO OFF CACHE BOOL "Use TAO in inverse solvers." FORCE)
#    message(STATUS  "TAO not found. Inverse solvers using the TAO library will not be built.")
#  endif()

  # Use option values to set compiler and linker flags
  set (Glint2_EXTERNAL_LIBS "")

#  # optional
#  if (Glint2_USE_PROJ4)
    include_directories (${PROJ4_INCLUDES})
    list (APPEND Glint2_EXTERNAL_LIBS ${PROJ4_LIBRARIES})
#  endif()

#  if (Glint2_USE_PNETCDF)
#    include_directories (${PNETCDF_INCLUDES})
#    list (APPEND Glint2_EXTERNAL_LIBS ${PNETCDF_LIBRARIES})
#  endif()

#  if (Glint2_USE_TAO)
#    include_directories (${TAO_INCLUDE_DIRS})
#    list (APPEND Glint2_EXTERNAL_LIBS ${TAO_LIBRARIES})
#  endif()

endmacro()

macro(glint2_set_dependencies)

  # Set include and library directories for *required* libraries.
  include_directories (
    ${PETSC_INCLUDES}
    ${FFTW_INCLUDE_DIRS}
    ${FFTW_INCLUDES}
    ${GSL_INCLUDES}
    ${UDUNITS2_INCLUDES}
    ${NETCDF_INCLUDES}
    ${MPI_C_INCLUDE_PATH})

  list (APPEND Glint2_EXTERNAL_LIBS
    ${PETSC_LIBRARIES}
    ${UDUNITS2_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${GSL_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${MPI_C_LIBRARIES})

  # Hide distracting CMake variables
  mark_as_advanced(file_cmd MPI_LIBRARY MPI_EXTRA_LIBRARY
    CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT
    MAKE_EXECUTABLE TAO_DIR TAO_INCLUDE_DIRS NETCDF_PAR_H)

endmacro()

macro(glint2_set_dependencies)

include_directories(
	${PROJECT_SOURCE_DIR}/slib
	${Boost_INCLUDE_DIRS}
    ${NETCDF_INCLUDES} ${NETCDF_INCLUDES_CXX} ${NETCDF_INCLUDES_F77}
    ${BLITZ++_INCLUDE_DIR}
    ${GMP_INCLUDE_DIR}
    ${CGAL_INCLUDE_DIR}
	# These should be optional
	${PISM_INCLUDE_DIRS}
	${GALAHAD_INCLUDE_DIRS})

list (APPEND Glint2_EXTERNAL_LIBS
	${Boost_LIBRARIES}
	${NETCDF_LIBRARIES}
	${BLITZ++_LIBRARY}
	${GMP_LIBRARY}
	${CGAL_LIBRARY}
	# These should be optional
	${PISM_LIBRARIES}
	${GALAHAD_LIBRARIES})

endmacro()
