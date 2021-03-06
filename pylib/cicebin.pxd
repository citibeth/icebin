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

    cdef cppclass RegridMatrices:
        RegridMatrices(RegridMatrices) except +

    cdef cppclass GCMRegridder:
        GCMRegridder() except +

        unsigned int nhc() except +
        unsigned long nA() except +
        unsigned long nE() except +

        void ncio(cibmisc.NcIO &, string vname) except +
        IceRegridder *ice_regridder(string name) except +

    cdef cppclass GCMRegridder_Standard(GCMRegridder):
        GCMRegridder_Standard() except +

cdef extern from "icebin_cython.hpp" namespace "icebin::cython":
    cdef void read_fgrid(
        cibmisc.unique_ptr[Grid] &fgridA,
        string &gridA_fname,
        string &gridA_vname) except +

    cdef cibmisc.shared_ptr[GCMRegridder] new_GCMRegridder_Standard(
        Grid &fgridA,
        vector[double] &hcdefs,
        bool correctA) except +

    cdef cibmisc.shared_ptr[GCMRegridder] new_GCMRegridder_WrapE(
        string &global_ecO,
        cibmisc.shared_ptr[GCMRegridder] &gcmO) except +

    cdef void GCMRegridder_WrapE_set_focean(
        GCMRegridder *_gcmA,
        PyObject *foceanAOp_py,
        PyObject *foceanAOm_py) except +

    cdef object GCMRegridder_wA(
        GCMRegridder *gcm_regridder,
        string &sheet_name,
        bool native,
        double fill) except +

    cdef void GCMRegridder_add_sheet(
        GCMRegridder *cself, Grid &fgridA,
        string &name,
        string &gridI_fname, string &gridI_vname,
        string &exgrid_fname, string &exgrid_vname,
        string &sinterp_style) except +

    cdef cibmisc.linear_Weighted *RegridMatrices_matrix(
        RegridMatrices *self, string spec_name) except +

    cdef Hntr_regrid(Hntr *hntr, object WTA_py, object A_py, bool mean_polar) except +

    cdef RegridMatrices *new_regrid_matrices(GCMRegridder *gcm, string &sheet_name, PyObject *elevmaskI_py,
        bool scale, bool correctA,
        double sigma_x, double sigma_y, double sigma_z, bool conserve) except +

    cdef object read_elevmask(string &xfname) except +

cdef extern from "icebin/GridSpec.hpp" namespace "icebin":
    pass

    cdef cppclass HntrSpec:
        const int im
        const int jm
        const double offi
        const double dlat

        HntrSpec(int, int, double, double)
        int size()
#        double dxyp(int)

cdef extern from "icebin/modele/hntr.hpp" namespace "icebin::modele":
    pass

    cdef cppclass Hntr:
        pass

        Hntr(double yp17, HntrSpec &A, HntrSpec &B, double DATMIS)

