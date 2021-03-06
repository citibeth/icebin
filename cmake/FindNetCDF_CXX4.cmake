# Input Variables
#    NETCDF_CXX4_ROOT
# Produces:
#    NETCDF_CXX4_LIBRARY
#    NETCDF_CXX4_INCLUDE_DIR


FIND_PATH(NETCDF_CXX4_INCLUDE_DIR netcdf
	HINTS ${NETCDF_CXX4_ROOT}/include)

FIND_LIBRARY(NETCDF_CXX4_LIBRARY NAMES netcdf_c++4 netcdf-cxx4
	HINTS ${NETCDF_CXX4_ROOT}/lib ${NETCDF_CXX4_ROOT}/lib64)


IF (NETCDF_CXX4_INCLUDE_DIR AND NETCDF_CXX4_LIBRARY)
   SET(NETCDF_CXX4_FOUND TRUE)
ENDIF (NETCDF_CXX4_INCLUDE_DIR AND NETCDF_CXX4_LIBRARY)

IF (NETCDF_CXX4_FOUND)
   IF (NOT NETCDF_CXX4_FIND_QUIETLY)
      MESSAGE(STATUS "Found NETCDF_CXX4_LIBRARY: ${NETCDF_CXX4_LIBRARY}")
   ENDIF (NOT NETCDF_CXX4_FIND_QUIETLY)
ELSE (NETCDF_CXX4_FOUND)
   IF (NETCDF_CXX4_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find NETCDF_CXX4")
   ENDIF (NETCDF_CXX4_FIND_REQUIRED)
ENDIF (NETCDF_CXX4_FOUND)
