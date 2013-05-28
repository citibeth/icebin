# this is example-file: src/bin/Makefile.am

include ../commonmk.am

AM_LDFLAGS += @BOOST_FILESYSTEM_LIBS@

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

#-L$(top_builddir)/slib -lglint2

AM_LDFLAGS += -L$(DESTDIR)$(libdir) -lglint2 \
	$(PYTHON_LDFLAGS) $(NUMPY_LDFLAGS)

# the program to build (the names of the final binaries)
bin_PROGRAMS = \
	hires \
	searise_a searise_g searise50 searise100 \
	greenland_4x5 greenland_2x2_5 \
	ga_2x2_5 \
	overlap \
	apitest blitztest smulttest

hires_SOURCES = hires.cpp

searise_g_SOURCES = searise_g.cpp
searise_a_SOURCES = searise_a.cpp

searise50_SOURCES = searise50.cpp
searise100_SOURCES = searise100.cpp

greenland_4x5_SOURCES = greenland_4x5.cpp
greenland_2x2_5_SOURCES = greenland_2x2_5.cpp

ga_2x2_5_SOURCES = greenland_2x2_5.cpp

overlap_SOURCES = overlap.cpp

apitest_SOURCES = apitest.cpp

test_grid2_SOURCES = test_grid2.cpp
test_grid1_SOURCES = test_grid1.cpp

blitztest_SOURCES = blitztest.cpp

smulttest_SOURCES = smulttest.cpp

###############################################################################
