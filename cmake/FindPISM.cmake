# - Try to find Pism
# Input Variables
#    PISM_DIR - Root of PISM install tree
# Once done this will define
#  PISM_FOUND - System has Pism
#  PISM_INCLUDE_DIRS - The Pism include directories
#  PISM_LIBRARIES - The libraries needed to use Pism
##  PISM_DEFINITIONS - Compiler switches required for using Pism

#pism_find_prerequisites()


find_path(PISM_INCLUDE_DIR pism/base/iceModel.hh
          HINTS ${PISM_DIR}/include)

set(PISM_INCLUDE_DIRS ${PISM_INCLUDE_DIR} ${PISM_INCLUDE_DIR}/pism)

find_library(PISM_LIBRARY NAMES pism
    HINTS ${PISM_DIR}/lib)
find_library(PISMICEBIN_LIBRARY NAMES pismicebin
    HINTS ${PISM_DIR}/lib)

set(PISM_LIBRARIES ${PISM_LIBRARY} ${PISMICEBIN_LIBRARY})
