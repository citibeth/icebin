# this is example-file: src/bin/Makefile.am

include ../commonmk.am

AM_LDFLAGS += @BOOST_FILESYSTEM_LIBS@

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

AM_LDFLAGS += $(top_builddir)/slib/libglint2.la

# the program to build (the names of the final binaries)
bin_PROGRAMS = #searise_ll_overlap searise_cesm_overlap searise_geos_overlap xy_overlap pytypes read_eqp

## list of sources for the 'helloWorld' binary
#searise_ll_overlap_SOURCES = searise_ll_overlap.cpp
#
## the additional libraries needed to link helloWorld
##searise_ll_overlap_LDADD = 	$(top_builddir)/slib/libglint2.la $(AM_LDFLAGS)
#
#searise_cesm_overlap_SOURCES = searise_cesm_overlap.cpp
#
#searise_geos_overlap_SOURCES = searise_geos_overlap.cpp
#
#xy_overlap_SOURCES = xy_overlap.cpp
#
##hctest_SOURCES = hctest.cpp
#
#pytypes_CXXFLAGS = $(AM_CXXFLAGS) $(PYTHON_CFLAGS) $(NUMPY_CFLAGS)
#pytypes_SOURCES = pytypes.cpp
#
#read_eqp_SOURCES = read_eqp.f90

###############################################################################
