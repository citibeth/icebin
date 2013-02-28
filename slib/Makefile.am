# this is example-file: src/glint2/Makefile.am

include ../commonmk.am

###############################################################################
# THE LIBRARIES TO BUILD
###############################################################################

# the library names to build (note we are building static libs only)
lib_LTLIBRARIES = libglint2.la

# the sources to add to the library and to add to the source distribution
libglint2_la_SOURCES = \
	giss/SparseMatrix.cpp \
	giss/hsl_zd11_x.cpp giss/ncutil.cpp giss/IndexTranslator.cpp \
	giss/geodesy.cpp \
	giss/Proj2.cpp \
	giss/f90blitz_f.f90 \
	glint2/Grid.cpp \
	glint2/Grid_LonLat.cpp \
	glint2/Grid_XY.cpp \
	glint2/read_grid.cpp \
	glint2/modele/grids_ll.cpp \
	glint2/gridutil.cpp \
	glint2/HeightClassifier.cpp \
	glint2/matrix_ops.cpp \
	glint2/IceSheet.cpp \
	glint2/IceSheet_L0.cpp \
	glint2/MatrixMaker.cpp \
	glint2/read_matrixmaker.cpp \
	glint2/ExchangeGrid.cpp \
	glint2/clippers.cpp

# # ------------------------------------------------------
giss/f90blitz_f.f90 ../include/giss/f90blitz_f.f90 : giss/f90blitz_f.py
	python giss/f90blitz_f.py giss/f90blitz_f.f90 ../include/giss/f90blitz_f.h

-include f90_deps.mk
FORTRAN_COMMAND_LINE = $(FCCOMPILE) -c -o

deps : f90_deps.mk
.PHONY : all deps

f90_deps.mk : ../f90_mod_deps.py $(FORTRAN_SOURCES)
	python ../f90_mod_deps.py $(FORTRAN_SOURCES) >f90_deps.mk

FORTRAN_SOURCES = giss/f90blitz_f.f90
libglint2_la_SOURCES += $(FORTRAN_SOURCES)

# Fortran headers (.mod files) are installed from this directory.
# For C/C++ header files, see ../include

#libglint2_la_HEADERS = $(filter-out eqp_f.mod,$(FORTRAN_SOURCES:.f90=.mod))
# We should really use a Python script to figure out which module files are created.
# See ../f90_mod_deps.py
libglint2_la_HEADERS = f90blitz.mod
libglint2_ladir = $(includedir)
# # ------------------------------------------------------


# PY_SOURCES = pyutil.cpp Grid_py.cpp Rasterizer_py.cpp Glint2_py.cpp glint2_module_py.cpp


# libglint2_la_SOURCES = $(shell find -E . -regex '.*\.(cpp|c|f90|F90|f|F)')
						
###############################################################################
