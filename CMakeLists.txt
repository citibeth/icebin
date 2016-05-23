cmake_minimum_required(VERSION 3.1)
project (GLINT2)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include (setup_rpath)
include_directories($ENV{CMAKE_TRANSITIVE_INCLUDE_PATH})

# ------- Process Configure Options

# -------- USE_PISM
if (NOT DEFINED USE_PISM)
	set(USE_PISM NO)
endif()
if(USE_PISM)
	add_definitions(-DUSE_PISM)
endif()

# ------- USE_PYTHON
if (NOT DEFINED USE_PYTHON)
	set(USE_PYTHON NO)
endif()

# -------- USE_EVERYTRACE
if (NOT DEFINED USE_EVERYTRACE)
	set(USE_EVERYTRACE YES)
endif()
if (USE_EVERYTRACE)
    add_definitions(-DUSE_EVERYTRACE)
endif()

# -------- Find Dependencies

# Used directly by Glint2
find_package(GALAHAD REQUIRED)		# This should be optional
find_package(Blitz++ REQUIRED)
find_package(CGAL REQUIRED)
find_package(MPFR REQUIRED)		# Used by CGAL
find_package(MPI REQUIRED)

find_package(NetCDF_C REQUIRED)
find_package(NetCDF_CXX REQUIRED)
find_package(NetCDF_Fortran REQUIRED)

find_package(Everytrace)

find_package(PROJ4)
find_package(Boost COMPONENTS filesystem system date_time thread REQUIRED)
find_package(UDUNITS2 REQUIRED)

# Used by dependencies
find_package(GMP REQUIRED)
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)

find_package(OpenMP REQUIRED)  # Used by GALAHAD, but not in code we run.  We are an MPI shop.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")


if (${USE_PISM})
	message('------------ Using PISM')
	find_package(PISM)
	find_package(PETSc)
endif()

if (${USE_PYTHON})
	message('------------ Using Python')
	find_package(Python3 REQUIRED)
	find_package(Numpy REQUIRED)
endif()



# Set Pism_EXTERNAL_LIBS and include directories.
include_directories (
	${EVERYTRACE_INCLUDE_DIR}
	${PROJECT_SOURCE_DIR}/slib
	${Boost_INCLUDE_DIRS}
    ${NETCDF_C_INCLUDE_DIR} ${NETCDF_CXX_INCLUDE_DIR} ${NETCDF_FORTRAN_INCLUDE_DIR}
    ${BLITZ++_INCLUDE_DIR}
    ${GMP_INCLUDE_DIR}
    ${CGAL_INCLUDE_DIR}
    ${FFTW_INCLUDE_DIRS}
    ${FFTW_INCLUDES}
    ${GSL_INCLUDES}
    ${UDUNITS2_INCLUDES}
    ${MPI_C_INCLUDE_PATH}
    ${PROJ4_INCLUDES}
	${GALAHAD_INCLUDE_DIRS})     # This should be optional

if (${USE_PISM})
	message(PISM_INCLUDE_DIRS ${PISM_INCLUDE_DIRS})
	include_directories(
		${PISM_INCLUDE_DIRS}
		${PETSC_INCLUDES}
	)
endif()

list(APPEND EXTERNAL_LIBS
	${EVERYTRACE_LIBRARY}
	${Boost_LIBRARIES}
    ${UDUNITS2_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${GSL_LIBRARIES}

    ${NETCDF_C_LIBRARY} ${NETCDF_CXX_LIBRARY} ${NETCDF_FORTRAN_LIBRARY}
    ${MPI_C_LIBRARIES}
	${GALAHAD_LIBRARIES}
	${PROJ4_LIBRARIES}
	${GMP_LIBRARY}
	${CGAL_LIBRARY}
	${MPFR_LIBRARIES}
)

if (${USE_PISM})
	list(APPEND EXTERNAL_LIBS
	    ${PETSC_LIBRARIES}
		${PISM_LIBRARIES}
	)
endif()


message(EXTERNAL_LIBS "${EXTERNAL_LIBS}")

# -------- Process subdirectories of the build
add_subdirectory(slib)
add_subdirectory(sbin)

if(${USE_PYTHON})
	add_subdirectory(pyext)
endif()

