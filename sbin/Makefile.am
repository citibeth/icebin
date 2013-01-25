# this is example-file: src/bin/Makefile.am

include ../commonmk.am

AM_LDFLAGS += @BOOST_FILESYSTEM_LIBS@

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

AM_LDFLAGS += $(top_builddir)/slib/libglint2.la

# the program to build (the names of the final binaries)
bin_PROGRAMS = searise_grid overlap test_grid1 test_grid2 greenland_4x5

searise_grid_SOURCES = searise_grid.cpp

overlap_SOURCES = overlap.cpp

greenland_4x5_SOURCES = greenland_4x5.cpp

test_grid2_SOURCES = test_grid2.cpp
test_grid1_SOURCES = test_grid1.cpp

###############################################################################
