# - Try to find Pism
# Input Variables
#    SPATIALINDEX_DIR - Root of SPATIALINDEX install tree
# Once done this will define
#  SPATIALINDEX_FOUND - System has Pism
#  SPATIALINDEX_INCLUDE_DIRS - The Pism include directories
#  SPATIALINDEX_LIBRARIES - The libraries needed to use Pism
##  SPATIALINDEX_DEFINITIONS - Compiler switches required for using Pism

#pism_find_prerequisites()


find_path(SPATIALINDEX_INCLUDE_DIR spatialindex/SpatialIndex.h
          HINTS ${SPATIALINDEX_DIR}/include)

set(SPATIALINDEX_INCLUDE_DIRS ${SPATIALINDEX_INCLUDE_DIR} ${SPATIALINDEX_INCLUDE_DIR}/spatialindex)

find_library(SPATIALINDEX_LIBRARY NAMES spatialindex
    HINTS ${SPATIALINDEX_DIR}/lib)
find_library(SPATIALINDEX_C_LIBRARY NAMES spatialindex_c
    HINTS ${SPATIALINDEX_DIR}/lib)

set(SPATIALINDEX_LIBRARIES ${SPATIALINDEX_LIBRARY} ${SPATIALINDEX_C_LIBRARY})
message("Found SPATIALINDEX_LIBRARIES ${SPATIALINDEX_LIBRARIES}")
