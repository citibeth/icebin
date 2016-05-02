cmake_minimum_required(VERSION 3.1)
project (GLINT2)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include (setup_rpath)
include_directories($ENV{CMAKE_TRANSITIVE_INCLUDE_PATH})

# ------- Process Configure Options

if (NOT DEFINED USE_EVERYTRACE)
    set(USE_EVERYTRACE NO)
endif()
if (USE_EVERYTRACE)
    find_package(Everytrace)
endif()


if (NOT DEFINED BUILD_PYTHON)
    set(BUILD_PYTHON NO)
endif()
if (BUILD_PYTHON)
    find_package(Python REQUIRED)
    find_package(Numpy REQUIRED)
endif()



if (NOT DEFINED BUILD_GRIDGEN)
    set(BUILD_GRIDGEN YES)
endif()
if (BUILD_GRIDGEN)
    find_package(CGAL REQUIRED)
    find_package(GMP REQUIRED)
    find_package(MPFR REQUIRED)

    list(APPEND EXTERNAL_LIBS ${CGAL_LIBRARY} ${GMP_LIBRARY} ${MPFR_LIBRARIES})
    list(APPEND Boost_COMPONENTS thread system program_options)
    include_directories(${CGAL_INCLUDE_DIR} ${GMP_INCLUDE_DIR} ${MPFR_INCLUDES})
    add_definitions(-DBUILD_GRIDGEN)
endif()

if (NOT DEFINED BUILD_COUPLER)
    set(BUILD_COUPLER YES)
endif()
if (BUILD_COUPLER)
    find_package(MPI REQUIRED)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    list(APPEND EXTERNAL_LIBS ${MPI_CXX_LIBRARIES})

    if (NOT DEFINED USE_PISM)
        set(USE_PISM NO)
    endif()
    if(USE_PISM)
        find_package(PISM)
        find_package(PETSC)
        include_directories(
            ${PISM_INCLUDE_DIRS}
            ${PETSC_INCLUDES}
        )
        list(APPEND EXTERNAL_LIBS ${PETSC_LIBRARIES} ${PISM_LIBRARIES})

        add_definitions(-DUSE_PISM)
    endif()
    add_definitions(-DBUILD_COUPLER)
endif()

find_package(Boost REQUIRED COMPONENTS ${Boost_COMPONENTS} filesystem date_time)
#find_package(Boost REQUIRED thread system log filesystem date_time)
include_directories(${Boost_INCLUDE_DIRS})
list(APPEND EXTERNAL_LIBS ${Boost_LIBRARIES})

# -------- Find Dependencies

# --- Blitz
find_package(Blitz REQUIRED)
add_definitions(-DUSE_BLITZ)
include_directories(${BLITZ_INCLUDE_DIR})
list(APPEND EXTERNAL_LIBS ${BLITZ_LIBRARY})

# --- NetCDF-CXX4
find_package(NetCDF_CXX4 REQUIRED)
add_definitions(-DUSE_NETCDF)
include_directories(${NETCDF_CXX4_INCLUDE_DIR})
list(APPEND EXTERNAL_LIBS ${NETCDF_CXX4_LIBRARY})


# --- IBMisc
find_package(Ibmisc REQUIRED)
list(APPEND EXTERNAL_LIBS ${IBMISC_LIBRARY})
include_directories(${IBMISC_INCLUDE_DIR})


# --- Proj.4
find_package(PROJ4 REQUIRED)
add_definitions(-DUSE_PROJ4)
include_directories(${PROJ4_INCLUDES})
list(APPEND EXTERNAL_LIBS ${PROJ4_LIBRARIES})

# --- Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

#find_package(GALAHAD REQUIRED)     # This should be optional

#find_package(CGAL REQUIRED)
#find_package(MPFR REQUIRED)        # Used by CGAL


#find_package(PROJ4)
#find_package(Boost COMPONENTS filesystem system date_time thread REQUIRED)
#find_package(UDUNITS2 REQUIRED)

## Used by dependencies
#find_package(GMP REQUIRED)
#find_package(GSL REQUIRED)
#find_package(FFTW REQUIRED)

#find_package(OpenMP REQUIRED)  # Used by GALAHAD, but not in code we run.  We are an MPI shop.
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")



## Set Pism_EXTERNAL_LIBS and include directories.
#include_directories (
#   ${EVERYTRACE_INCLUDE_DIR}
#   ${PROJECT_SOURCE_DIR}/slib
#   ${Boost_INCLUDE_DIRS}
#    ${NETCDF_INCLUDES} ${NETCDF_INCLUDES_CXX} ${NETCDF_INCLUDES_F77}
#    ${BLITZ++_INCLUDE_DIR}
#    ${GMP_INCLUDE_DIR}
#    ${CGAL_INCLUDE_DIR}
#    ${FFTW_INCLUDE_DIRS}
#    ${FFTW_INCLUDES}
#    ${GSL_INCLUDES}
#    ${UDUNITS2_INCLUDES}
#    ${MPI_C_INCLUDE_PATH}
#
#   ${GALAHAD_INCLUDE_DIRS})     # This should be optional


#list(APPEND EXTERNAL_LIBS
#   ${EVERYTRACE_CF_LIBRARY}
#   ${Boost_LIBRARIES}
#    ${UDUNITS2_LIBRARIES}
#    ${FFTW_LIBRARIES}
#    ${GSL_LIBRARIES}
#
#    ${NETCDF_LIBRARIES_F90} ${NETCDF_LIBRARIES_CXX} ${NETCDF_LIBRARIES_C}
#    ${MPI_C_LIBRARIES}
#   ${GALAHAD_LIBRARIES}
#   ${PROJ4_LIBRARIES}
#   ${GMP_LIBRARY}
#   ${CGAL_LIBRARY}
#   ${MPFR_LIBRARIES}
#)
#
#if (${USE_PISM})
#   list(APPEND EXTERNAL_LIBS
#       ${PETSC_LIBRARIES}
#       ${PISM_LIBRARIES}
#   )
#endif()

# -------- Process subdirectories of the build
add_subdirectory(slib)
#add_subdirectory(sbin)
add_subdirectory(tests)
if (BUILD_GRIDGEN)
    add_subdirectory(grids)
endif()

if(${BUILD_PYTHON})
    add_subdirectory(pylib)
endif()

