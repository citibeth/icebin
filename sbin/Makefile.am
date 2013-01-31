# this is example-file: src/bin/Makefile.am

include ../commonmk.am

AM_LDFLAGS += @BOOST_FILESYSTEM_LIBS@

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

AM_LDFLAGS += $(top_builddir)/slib/libglint2.la \
	$(PYTHON_LDFLAGS) $(NUMPY_LDFLAGS)

# the program to build (the names of the final binaries)
bin_PROGRAMS = \
	hires \
	searise searise50 \
	greenland_4x5 greenland_2x2_5 \
	overlap

hires_SOURCES = hires.cpp

searise_SOURCES = searise.cpp
searise50_SOURCES = searise50.cpp

greenland_4x5_SOURCES = greenland_4x5.cpp
greenland_2x2_5_SOURCES = greenland_2x2_5.cpp

overlap_SOURCES = overlap.cpp

test_grid2_SOURCES = test_grid2.cpp
test_grid1_SOURCES = test_grid1.cpp

###############################################################################
