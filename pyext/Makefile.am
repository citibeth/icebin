
include ../commonmk.am

# PYTHON_PREFIX := $(shell $(PYTHON)-config --prefix)

AM_CXXFLAGS += $(PYTHON_CFLAGS) $(NUMPY_CFLAGS)
pyexec_LTLIBRARIES = glint2.la
glint2_la_SOURCES = #pyutil.cpp Grid_py.cpp Rasterizer_py.cpp Glint2_py.cpp glint2_module_py.cpp
glint2_la_LDFLAGS = -module $(PYTHON_LDFLAGS) $(NUMPY_LDFLAGS) $(AM_LDFLAGS) $(top_srcdir)/slib/libglint2.la

