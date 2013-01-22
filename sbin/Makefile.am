# this is example-file: src/bin/Makefile.am

include ../commonmk.am

AM_LDFLAGS += @BOOST_FILESYSTEM_LIBS@

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

AM_LDFLAGS += $(top_builddir)/slib/libglint2.la

# the program to build (the names of the final binaries)
bin_PROGRAMS = searise_grid

searise_grid_SOURCES = searise_grid.cpp

###############################################################################
