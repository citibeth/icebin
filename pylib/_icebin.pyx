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

# On reference counting...
# https://groups.google.com/forum/#!topic/cython-users/pZFurj3rUyg

from cpython.object cimport *       # PyObject
cimport cicebin
cimport cibmisc     # C++ stuff in ibmisc
cimport ibmisc      # Cython stuff in ibmisc
import numpy as np
cimport numpy as np
np.import_array()
import scipy.sparse
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp cimport bool
import functools
import operator

cdef class IceRegridder:
    pass


cdef split_shape(ashape, alen):
    # See if we can find some set of dimensions matching ilen
    cumlen = np.cumprod(tuple(reversed(ashape)))
    try:
        icut = len(ashape) - 1 - next(i for i,x in enumerate(cumlen) if x==alen)
    except StopIteration:
        raise ValueError('Cannot find trailing dimension of {} in input of shape {}'.format(alen, ashape)) from None

    return ashape[0:icut], (
        functools.reduce(operator.mul, ashape[0:icut], 1),
        functools.reduce(operator.mul, ashape[icut:], 1))


cdef class WeightedSparse:
    """The result of RegridMatrices.matrix()"""
    cdef cicebin.CythonWeightedSparse *cself

    def __dealloc__(self):
        del self.cself

    def __call__(self):
        """Obtain the matrix and weight vectors as Python structures.
        returns: (wM, M, Mw)
            wM: np.array
                Weight vector ("area" of output grid cells)
            M: scipy.sparse.coo_matrix
               Regridding Matrix
            Mw: np.array
                Weight vector ("area" of input grid cells)
       """

        wM, (data,shape), Mw = cicebin.CythonWeightedSparse_to_tuple(self.cself)
        return wM, scipy.sparse.coo_matrix(data, shape), Mw

    def apply(self, A_s, double fill=np.nan):
        """Applies the regrid matrix to A_s.  Smoothes and conserves, if those
        options were specified in RegridMatrices.matrix().
        A_s: Either:
            - A single vector (1-D array) to be transformed.
            - A 2-D array of row vectors to be transformed.
        fill:
            Un-set indices in output array will get this value."""

        # Number of elements in sparse in put vector
        _,alen = cicebin.CythonWeightedSparse_sparse_extent(self.cself)

        leading_shape, new_shape = split_shape(A_s.shape, alen)
        A_s = A_s.reshape(new_shape)
        B_s = cicebin.CythonWeightedSparse_apply(self.cself, <PyObject *>A_s, fill)

        B_s = B_s.reshape( leading_shape + (B_s.shape[1],) )

        return B_s

cdef class RegridMatrices:
    cdef cicebin.RegridMatrices *cself

    def __dealloc__(self):
        del self.cself

    def matrix(self, str spec_name, bool scale=True, bool correctA=True,
        sigma=(0,0,0), conserve=True):
        """Compute a regrid matrix.
        spec_name:
            Type of regrid matrix to obtain.  Choice are:
            'EvI', 'AvI', 'IvA', 'IvE', 'EvA', 'AvE'
        scale: Produce scaled matrix?
            true  --> [kg m-2]
            false --> [kg]
        correctA: Correct for projection on A and E side of matrices?
        conserve: Apply conservation correction
            (if needed; eg, if sigma is set, and smoothing is supported
            for the output grid).
        returns: WeightedSparse
        """
        cdef cicebin.CythonWeightedSparse *crm
        crm = cicebin.RegridMatrices_matrix(
            self.cself, spec_name.encode(), scale, correctA,
            sigma[0], sigma[1], sigma[2], conserve)
        ret = WeightedSparse()
        ret.cself = crm
        return ret

cdef class GCMRegridder:
    cdef cicebin.GCMRegridder *cself

    def __dealloc__(self):
        del self.cself

    def __init__(self, *args):
        cdef ibmisc.NcIO ncio

        # Create a brand new GCMRegridder
        if len(args) == 4:
            gridA_fname, gridA_vname, hcdefs, correctA = args
            self.cself = cicebin.new_GCMRegridder_Standard(
                gridA_fname.encode(), gridA_vname.encode(), hcdefs, correctA)
        elif len(args) == 1:
            self.cself = new cicebin.GCMRegridder_Standard()
            (regridder_fname,) = args

            # Load an existing GCMRegridder from disk
            ncio = ibmisc.NcIO(regridder_fname, 'read')
            self.ncio(ncio, str('m'))
            ncio.close()
        else:
            raise ValueError('Invalid arguments: {}'.format(args))

    def to_modele(self):
        cdef cicebin.GCMRegridder *gcm
        gcm = cicebin.new_GCMRegridder_ModelE(self.cself)
        if not gcm:
            raise RuntimeError('IceBin must be built with USE_MODELE in order to use ModelE features')
        self.cself = gcm

    def ncio(self, ibmisc.NcIO ncio, vname):
        self.cself.ncio(deref(ncio.cself), vname.encode())

    def add_sheet(self, name,
        gridI_fname, gridI_vname,
        exgrid_fname, exgrid_vname,
        interp_style,
        elevI, maskI):

        elevI = elevI.reshape(-1)
        maskI = maskI.reshape(-1)
        elevI[maskI != 0] = np.nan  # For IceBin, isnan(elevI) means it's masked out.
        cicebin.GCMRegridder_add_sheet(self.cself,
            name.encode(),
            gridI_fname.encode(), gridI_vname.encode(),
            exgrid_fname.encode(), exgrid_vname.encode(),
            interp_style.encode(),
            <PyObject *>elevI)   # Borrowed references

    def set_elevI(self, sheet_name, elevI):
        elevI = elevI.reshape(-1)
        cicebin.GCMRegridder_set_elevI(self.cself,
            sheet_name.encode(), <PyObject *>elevI)    # Borrowed references

    def regrid_matrices(self, str sheet_name):
        cdef cicebin.RegridMatrices *crm = new cicebin.RegridMatrices(
            self.cself.regrid_matrices(sheet_name.encode()))

#        cdef cicebin.RegridMatrices *crm = new cicebin.RegridMatrices(
#            self.cself.regrid_matrices(sheet_name.encode()))
        rm = RegridMatrices()
        rm.cself = crm
        return rm

def coo_multiply(M, xx, double fill=np.nan, ignore_nan=False):
    """M:
        SciPy sparse matrix"""
    xx = xx.reshape(-1)
    yy = np.zeros(M._shape[0])
    yy[:] = fill

    cicebin.coo_matvec(<PyObject *>yy, <PyObject *>xx, ignore_nan,
        M._shape[0], M._shape[1],
        <PyObject *>M.row, <PyObject *>M.col, <PyObject *>M.data)

    return yy

# ============================================================

cdef class HntrGrid:
    cdef cicebin.HntrGrid *cself;

    def __dealloc__(HntrGrid self):
        del self.cself

    def __init__(HntrGrid self, int im, int jm, float offi, float dlat):
        self.cself = new cicebin.HntrGrid(im, jm, offi, dlat)

    @property
    def im(self):
        return self.cself.im

    @property
    def jm(self):
        return self.cself.jm

    @property
    def size(self):
        return self.cself.size()

    @property
    def offi(self):
        return self.cself.offi

    @property
    def dlat(self):
        return self.cself.dlat

    @property
    def dxyp(self, int j):
        return self.cself.dxyp(j)


cdef class Hntr:
    cdef cicebin.Hntr *cself;

    def __dealloc__(self):
        del self.cself

    def __init__(self, HntrGrid Agrid, HntrGrid Bgrid, float DATMIS):
        self.cself = new cicebin.Hntr(Agrid.cself[0], Bgrid.cself[0], DATMIS)

    def regrid(self, WTA, A, bool mean_polar):
        cicebin.Hntr_regrid(self.cself, WTA, A, mean_polar)
