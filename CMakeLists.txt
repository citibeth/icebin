# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (Glint2)
enable_language(Fortran)
include ("CMake/PISM_CMake_macros.cmake")
include ("CMake/GLINT2_CMake_macros.cmake")
list (APPEND CMAKE_MODULE_PATH "${Glint2_SOURCE_DIR}/CMake")

# ==========================================================
## http://stackoverflow.com/questions/2210839/conditional-cxx-flags-using-cmake-based-on-compiler
## Also, special rounding required for CGAL (-frounding-math)
if (CMAKE_COMPILER_IS_GNUCXX)
  set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -frounding-math"
	CACHE STRING "g++ Compiler Flags for Debug Builds" FORCE)
#
##  set ( CMAKE_CXX_FLAGS_DEBUG "/DEBUG_FLAGS_GO_HERE")
##  set ( CMAKE_CXX_FLAGS_RELEASE  "/RELEASE_FLAGS_GO_HERE" )
endif ()

# See:
## http://cmake.3232098.n2.nabble.com/CMake-rebuilding-too-much-td5611087.html
#if(CMAKE_COMPILER_IS_GNUCXX AND NOT __CHECKED_CXX_FLAGS)
#  # this is the guard so we do this only once
#  set(__CHECKED_CXX_FLAGS TRUE CACHE INTERNAL "Whether we checked the CXX flags already")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -frounding-math "
#	CACHE STRING "g++ Compiler Flags for Debug Builds" FORCE)
##	Wall -Woverloaded-virtual -Wno-unknown-pragmas -Werror"
##      CACHE STRING "Flags used by the compiler during all build types." FORCE)
##  endif()
#endif() 



# ==========================================================
set(NETCDF_F77 "YES")
set(NETCDF_CXX "YES")
#pism_find_prerequisites()
glint2_find_prerequisites()

message(NETCDF ${NETCDF_LIBRARIES})


find_package(GALAHAD)
find_package(Blitz++)
find_package(GMP)
find_package(CGAL)
find_package(Boost COMPONENTS filesystem system date_time)

find_package(PISM)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_PISM"
	CACHE STRING "g++ Compiler Flags for Debug Builds" FORCE)



MESSAGE(PISM ${PISM_INCLUDE_DIRS})
MESSAGE(PISM ${PISM_LIBRARIES})


# ==========================================================
# See "Always full RPATH"
# http://www.cmake.org/Wiki/CMake_RPATH_handling#CMake_and_the_RPATH

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")
# ==========================================================

## Set rpath in the executable to point to the library
## http://cmake.3232098.n2.nabble.com/Install-an-executable-with-his-libraries-td7199530.html
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib") 



# Recurse into the "Hello" and "Demo" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (slib)
add_subdirectory (sbin)
