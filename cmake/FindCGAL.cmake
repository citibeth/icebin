# Input Variables
#    CGAL_ROOT
# Produces:
#    CGAL_LIBRARY
#    CGAL_INCLUDE_DIR


FIND_PATH(CGAL_INCLUDE_DIR CGAL/version.h
	HINTS ${CGAL_ROOT}/include)

FIND_LIBRARY(CGAL_LIBRARY NAMES CGAL
	HINTS ${CGAL_ROOT}/lib)

IF (CGAL_INCLUDE_DIR AND CGAL_LIBRARY)
   SET(CGAL_FOUND TRUE)
ENDIF (CGAL_INCLUDE_DIR AND CGAL_LIBRARY)

IF (CGAL_FOUND)
   IF (NOT CGAL_FIND_QUIETLY)
      MESSAGE(STATUS "Found CGAL_LIBRARY: ${CGAL_LIBRARY}")
   ENDIF (NOT CGAL_FIND_QUIETLY)
ELSE (CGAL_FOUND)
   IF (CGAL_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find CGAL")
   ENDIF (CGAL_FIND_REQUIRED)
ENDIF (CGAL_FOUND)
