# - Try to find Pism
# Input Variables
#    PISM_SRC
#    PISM_LIB
# Once done this will define
#  PISM_FOUND - System has Pism
#  PISM_INCLUDE_DIRS - The Pism include directories
#  PISM_LIBRARIES - The libraries needed to use Pism
##  PISM_DEFINITIONS - Compiler switches required for using Pism

#pism_find_prerequisites()


find_path(PISM_SRC base/iceModel.hh
          HINTS ${PISM_SRC})

set(PISM_INCLUDE_DIRS ${PISM_SRC})
#  ${PISM_SRC}/base
#  ${PISM_SRC}/base/stressbalance
#  ${PISM_SRC}/base/util
#  ${PISM_SRC}/base/util/io
#  ${PISM_SRC}/base/energy
#  ${PISM_SRC}/base/rheology
#  ${PISM_SRC}/base/basalstrength
#  ${PISM_SRC}/earth
#  ${PISM_SRC}/coupler
#  ${PISM_SRC}/coupler/atmosphere
#  ${PISM_SRC}/coupler/surface
#  ${PISM_SRC}/coupler/ocean
#  ${PISM_SRC}/coupler/util)


#     -DUSE_PISM @PETSC_CFLAGS@)

set(PISM_COMPONENTS base earth boundary stressbalance flowlaws util)
foreach(COMPONENT ${PISM_COMPONENTS})
    string(TOUPPER ${COMPONENT} UPPERCOMPONENT)
	find_library(PISM_${UPPERCOMPONENT}_LIBRARY pism${COMPONENT}
		PATH_SUFFIXES pism
		HINTS ${PISM_LIB})
	set(PISM_LIBRARIES ${PISM_LIBRARIES} ${PISM_${UPPERCOMPONENT}_LIBRARY})
endforeach()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PISM_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Pism  DEFAULT_MSG
                                  PISM_LIBRARIES PISM_INCLUDE_DIRS)

mark_as_advanced(PISM_INCLUDE_DIRS PISM_LIBRARIES )
