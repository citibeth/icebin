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
import warnings

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

    @property
    def shape(self):
        return self.cself.shape()

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

    def apply(self, A_s, double fill=np.nan, force_conservation=True):
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
        B_s = cicebin.CythonWeightedSparse_apply(self.cself, <PyObject *>A_s, fill, force_conservation)

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
    cdef cibmisc.shared_ptr[cicebin.GCMRegridder] cself
    cdef cibmisc.unique_ptr[cicebin.Grid] fgridA

    def __init__(self, *args):
        cdef ibmisc.NcIO ncio

        # Create a brand new GCMRegridder
        if len(args) == 4:
            gridA_fname, gridA_vname, hcdefs, correctA = args
            
            cicebin.read_fgrid(self.fgridA, gridA_fname.encode(), gridA_vname.encode())
            self.cself = cicebin.new_GCMRegridder_Standard(self.fgridA.get()[0], hcdefs, correctA)
        elif len(args) == 1:
            self.cself.reset(new cicebin.GCMRegridder_Standard())
            (regridder_fname,) = args

            # Load an existing GCMRegridder from disk
            ncio = ibmisc.NcIO(regridder_fname, 'read')
            self.ncio(ncio, str('m'))
            ncio.close()
        elif len(args) == 0:
            pass
            # Initialize blank
            # self.cself.reset(NULL)
        else:
            raise ValueError('Invalid arguments: {}'.format(args))

    def wA(self, sheet_name, snative, fill=0.):
        """Returns weights (as a vector) of overall grid."""
        if snative == 'native':
            native = True
        elif snative == 'proj':
            native = False
        else:
            raise ValueError("Invalid argument: snative must be 'native' or 'proj'")

        return cicebin.GCMRegridder_wA(self.cself.get(), sheet_name.encode(), native, fill)


    def to_modele(self, focean=None):
        """Returns a new GCMRegridder object, suitable for use with ModelE

        focean: (foceanAOp, foceanAOm) [OPTIONAL]
            Use this ocean in GCMRegridder_ModelE."""
        cdef cibmisc.shared_ptr[cicebin.GCMRegridder] gcm
        gcm = cicebin.new_GCMRegridder_ModelE(self.cself)
        if not gcm.get():
            raise RuntimeError('IceBin must be built with USE_MODELE in order to use ModelE features')

        ret = GCMRegridder()
        ret.cself = gcm

        if focean is not None:
            foceanAOp,foceanAOm = focean
            foceanAOp = foceanAOp.reshape(-1)
            foceanAOm = foceanAOm.reshape(-1)
            cicebin.GCMRegridder_ModelE_set_focean(ret.cself.get(), <PyObject *>foceanAOp, <PyObject *>foceanAOm)

        return ret

    def update_topo(self, topoO_fname, elevmask_sigmas, bool initial_timestamp, segments, primary_segment,
        fhc, underice, elevE,
        focean, flake, fgrnd, fgice, zatmo,
        foceanOm0):
        """Allows Python users access to GCMCoupler_Modele::update_topo().
        Starting from output of Gary's program (on the Ocean grid), this subroutine
        produces a ModelE TOPO file (as internal arrays) on the Atmosphere grid."""

        cicebin.update_topo(
            self.cself.get(), topoO_fname.encode(),
            <PyObject *>elevmask_sigmas, initial_timestamp, segments.encode(), primary_segment.encode(),
            <PyObject *>fhc, <PyObject *>underice, <PyObject *>elevE,
            <PyObject *>focean, <PyObject *>flake, <PyObject *>fgrnd, <PyObject *>fgice, <PyObject *>zatmo,
            <PyObject *>foceanOm0)

    def ncio(self, ibmisc.NcIO ncio, vname):
        self.cself.get().ncio(deref(ncio.cself), vname.encode())

    def add_sheet(self, name,
        gridI_fname, gridI_vname,
        exgrid_fname, exgrid_vname,
        interp_style):

        cicebin.GCMRegridder_add_sheet(self.cself.get(),
            self.fgridA.get()[0],
            name.encode(),
            gridI_fname.encode(), gridI_vname.encode(),
            exgrid_fname.encode(), exgrid_vname.encode(),
            interp_style.encode())

    def regrid_matrices(self, str sheet_name, elevI):
        elevI = elevI.reshape(-1)
        cdef cicebin.RegridMatrices *crm = \
            cicebin.new_regrid_matrices(self.cself.get(), sheet_name.encode(),
            <PyObject *>elevI)   # PyObject=Borrowed reference, object = owned reference
        rm = RegridMatrices()
        rm.cself = crm
        return rm

def coo_multiply(M, xx, double fill=np.nan, ignore_nan=False):
    """M:
        SciPy sparse matrix"""

    warnings.warn(
        "coo_multiply is deprecated; use WeightedSparse.apply() instead.",
        DeprecationWarning)

    xx = xx.reshape(-1)
    yy = np.zeros(M._shape[0])
    yy[:] = fill

    cicebin.coo_matvec(<PyObject *>yy, <PyObject *>xx, ignore_nan,
        M._shape[0], M._shape[1],
        <PyObject *>M.row, <PyObject *>M.col, <PyObject *>M.data)

    return yy

# ============================================================

cdef class HntrSpec:
    cdef cicebin.HntrSpec *cself;

    def __dealloc__(HntrSpec self):
        del self.cself

    def __init__(HntrSpec self, int im, int jm, float offi, float dlat):
        self.cself = new cicebin.HntrSpec(im, jm, offi, dlat)

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

#    @property
#    def dxyp(self, int j):
#        return self.cself.dxyp(j)


cdef class Hntr:
    cdef cicebin.Hntr *cself;

    def __dealloc__(self):
        del self.cself

    def __init__(self, double yp17, HntrSpec Bgrid, HntrSpec Agrid, float DATMIS):
        self.cself = new cicebin.Hntr(17.17, Bgrid.cself[0], Agrid.cself[0], DATMIS)

    def regrid(self, WTA, A, bool mean_polar):
        cicebin.Hntr_regrid(self.cself, WTA, A, mean_polar)

