# - Try to find BLITZ++
# Once done, this will define
#
# BLITZ++_FOUND - system has BLITZ++
# BLITZ++_INCLUDE_DIR - the BLITZ++ include directory
# BLITZ++_LIBRARY - lib to link to use BLITZ++

include(LibFindMacros)

# Use pkg-config to get hints about paths
#libfind_pkg_check_modules(BLITZ++_PKGCONF blitz)

# Include dir
find_path(BLITZ++_INCLUDE_DIR
  NAMES blitz/blitz.h
  HINTS ${BLITZ++_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(BLITZ++_LIBRARY
  NAMES blitz
  HINTS ${BLITZ++_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(BLITZ++_PROCESS_INCLUDE BLITZ++_INCLUDE_DIR )
set(BLITZ++_PROCESS_LIB BLITZ++_LIBRARY)
libfind_process(BLITZ++)
