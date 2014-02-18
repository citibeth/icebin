# - Try to find Galahad
# Input Variables
#    GALAHAD_ROOT
#    GALAHAD_ARCH
# Once done this will define
#  GALAHAD_FOUND - System has Galahad
#  GALAHAD_INCLUDE_DIRS - The Galahad include directories
#  GALAHAD_LIBRARIES - The libraries needed to use Galahad
##  GALAHAD_DEFINITIONS - Compiler switches required for using Galahad

find_path(GALAHAD_ROOT  modules/${GALAHAD_ARCH}/double/galahad_qpt_double.mod
          HINTS ${GALAHAD_ROOT})

message(GALAHAD_ROOT ${GALAHAD_ROOT})

set(GALAHAD_INCLUDE_DIRS ${GALAHAD_ROOT}/modules/${GALAHAD_ARCH}/double)

set(GALAHAD_LIB ${GALAHAD_ROOT}/objects/${GALAHAD_ARCH}/double)

#     -DUSE_GALAHAD @PETSC_CFLAGS@)

set(GALAHAD_COMPONENTS galahad galahad_hsl galahad_pardiso galahad_wsmp galahad_metis galahad_lapack galahad_blas)
set(GALAHAD_LIBRARIES gomp)		# Part of GCC
foreach(COMPONENT ${GALAHAD_COMPONENTS})
    string(TOUPPER ${COMPONENT} UPPERCOMPONENT)
	find_library(${UPPERCOMPONENT}_LIBRARY ${COMPONENT}
		HINTS ${GALAHAD_LIB})
	set(GALAHAD_LIBRARIES ${GALAHAD_LIBRARIES} ${${UPPERCOMPONENT}_LIBRARY})
endforeach()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GALAHAD_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Galahad  DEFAULT_MSG
                                  GALAHAD_LIBRARIES GALAHAD_INCLUDE_DIRS)

mark_as_advanced(GALAHAD_INCLUDE_DIRS GALAHAD_LIBRARIES )
