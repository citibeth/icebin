cmake_minimum_required(VERSION 3.1)
project (GLINT2)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include (setup_rpath)

# ------- Process Configure Options
if (NOT DEFINED USE_PISM)
	set(USE_PISM NO)
endif()
if(${USE_PISM})
	add_definitions(-DUSE_PISM)
endif()

if (NOT DEFINED USE_PYTHON)
	set(USE_PYTHON NO)
endif()

# -------- Find Dependencies

# Used directly by Glint2
find_package(GALAHAD REQUIRED)		# This should be optional
find_package(Blitz++ REQUIRED)
find_package(CGAL REQUIRED)
find_package(MPFR REQUIRED)		# Used by CGAL
find_package(MPI REQUIRED)

set(NETCDF_CXX "YES")
set(NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)

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
	find_package(PETSC)
endif()

if (${USE_PYTHON})
	find_package(Python REQUIRED)
	find_package(Numpy REQUIRED)
endif()



# Set Pism_EXTERNAL_LIBS and include directories.
include_directories (
	${PROJECT_SOURCE_DIR}/slib
	${Boost_INCLUDE_DIRS}
    ${NETCDF_INCLUDES} ${NETCDF_INCLUDES_CXX} ${NETCDF_INCLUDES_F77}
    ${BLITZ++_INCLUDE_DIR}
    ${GMP_INCLUDE_DIR}
    ${CGAL_INCLUDE_DIR}
    ${FFTW_INCLUDE_DIRS}
    ${FFTW_INCLUDES}
    ${GSL_INCLUDES}
    ${UDUNITS2_INCLUDES}
    ${MPI_C_INCLUDE_PATH}

	${GALAHAD_INCLUDE_DIRS})     # This should be optional

if (${USE_PISM})
	message(PISM_INCLUDE_DIRS ${PISM_INCLUDE_DIRS})
	include_directories(
		${PISM_INCLUDE_DIRS}
		${PETSC_INCLUDES}
	)
endif()

message('zyx' ${GALAHAD_LIBRARIES})

list(APPEND EXTERNAL_LIBS
	${Boost_LIBRARIES}
    ${UDUNITS2_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${GSL_LIBRARIES}

    ${NETCDF_LIBRARIES}
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

