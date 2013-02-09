
include ../commonmk.am

# PYTHON_PREFIX := $(shell $(PYTHON)-config --prefix)

# AM_CXXFLAGS += $(PYTHON_CFLAGS) $(NUMPY_CFLAGS)

pyexec_LTLIBRARIES = _glint2.la
_glint2_la_SOURCES = pyutil.cpp Grid_py.cpp MatrixMaker_py.cpp matrix_ops_py.cpp _glint2_module_py.cpp

_glint2_la_LDFLAGS = -module $(PYTHON_LDFLAGS) $(NUMPY_LDFLAGS) $(AM_LDFLAGS) $(top_srcdir)/slib/libglint2.la
