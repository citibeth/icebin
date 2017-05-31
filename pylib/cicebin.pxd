# IceBin: A Coupling Library for Ice Models and GCMs
# Copyright (c) 2013-2016 by Elizabeth Fischer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from cpython.object cimport *   # PyObject
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
cimport cibmisc

# This will make the C++ class def for Rectangle available..
cdef extern from "icebin/Grid.hpp" namespace "icebin":
    pass

    cdef cppclass Vertex:
        Vertex(double, double, long) except +
        long index
        const double x
        const double y


    cdef cppclass Cell:
        Cell(vector[Vertex *]) except +

    cdef cppclass GridMap[T]:
        cppclass iterator:
            T operator*()
            bint operator==(iterator)
            bint operator!=(iterator)
        iterator begin()
        iterator end()
        T *at(long)
        T *add_claim(T *)

    cdef cppclass Grid:
        GridMap[Vertex] vertices
        GridMap[Cell] cells

        Grid() except +


cdef extern from "icebin/GCMRegridder.hpp" namespace "icebin":
    pass

    cdef cppclass IceRegridder:
        pass

    cdef cppclass GCMRegridder:
        GCMRegridder() except +
        void ncio(cibmisc.NcIO &, string vname) except +
        IceRegridder *ice_regridder(string name) except +


    cdef cppclass RegridMatrices:
        RegridMatrices(IceRegridder *) except +


cdef extern from "icebin_cython.hpp" namespace "icebin::cython":
    cdef cppclass CythonWeightedSparse:
        pass

    cdef void GCMRegridder_init(
        GCMRegridder *self,
        string &gridA_fname,
        string &gridA_vname,
        vector[double] &hcdefs,
        bool correctA) except +

    cdef void GCMRegridder_add_sheet(
        GCMRegridder *cself,
        string &name,
        string &gridI_fname, string &gridI_vname,
        string &exgrid_fname, string &exgrid_vname,
        string &sinterp_style,
        PyObject *elevI_py) except +        # PyObject=Borrowed reference, object = owned reference

    cdef void GCMRegridder_set_elevI(
        GCMRegridder *cself,
        string &name,
        PyObject *elevI_py) except +        # PyObject=Borrowed reference, object = owned reference

    cdef CythonWeightedSparse *RegridMatrices_matrix(
        RegridMatrices *self, string spec_name, bool scale, bool correctA,
        double sigma_x, double sigma_y, double sigma_z, bool conserve) except +

    cdef object CythonWeightedSparse_dense_extent(CythonWeightedSparse *self) except +
    cdef object CythonWeightedSparse_sparse_extent(CythonWeightedSparse *self) except +
    cdef object CythonWeightedSparse_to_tuple(CythonWeightedSparse *self) except +

    cdef object CythonWeightedSparse_apply(CythonWeightedSparse *BvA, PyObject *A, double fill) except +

    cdef void coo_matvec(PyObject *yy_py, PyObject *xx_py, bool ignore_nan,
        int M_nrow, int M_ncol, PyObject *M_row_py, PyObject *M_col_py, PyObject *M_data_py) except +
