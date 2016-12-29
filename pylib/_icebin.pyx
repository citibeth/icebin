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

cdef class IceRegridder:
    pass

cdef class RegridMatrices:
    cdef cicebin.RegridMatrices *cself

    def __dealloc__(self):
        del self.cself

    def regrid(self, str spec_name, scale=True, correctA=True):
        """Compute a regrid matrix.
        spec_name:
            Type of regrid matrix to obtain.  Choice are:
            'EvI', 'AvI', 'IvA', 'IvE', 'EvA', 'AvE'
        scale: Produce scaled matrix?
            true  --> [kg m-2]
            false --> [kg]
        correctA: Correct for projection on A and E side of matrices?
        returns: (M, weights)
            M: scipy.sparse.coo_matrix
                Unscaled regridding matrix (i.e. produces [kg] not [kg m-2])
            weight: np.array
                Weight vector, defined by:
                    M(scaled=True) = diag(1/weight) * M(scaled=False)
                    M(scaled=False) = diag(weight) * M(scaled=True)
        """
        (data,shape), weight = cicebin.RegridMatrices_regrid(self.cself, spec_name.encode(), scale, correctA)
        # scipy.sparse.coo_matrix((data1, (rows1, cols1)), shape=(nrow1, ncol1))
        return scipy.sparse.coo_matrix(data, shape), weight

cdef class GCMRegridder:
    cdef cicebin.GCMRegridder cself

    def __init__(self, *args):
        cdef ibmisc.NcIO ncio

        try:
            gridA_fname, gridA_vname, hcdefs, correctA = args
            cicebin.GCMRegridder_init(&self.cself, gridA_fname.encode(), gridA_vname.encode(), hcdefs, correctA)
            return
        except:
            pass

        (regridder_fname,) = args
        ncio = ibmisc.NcIO(regridder_fname, 'read')
        self.ncio(ncio, str('m'))
        ncio.close()

    def ncio(self, ibmisc.NcIO ncio, vname):
        self.cself.ncio(deref(ncio.cself), vname.encode())

    def add_sheet(self, name,
        gridI_fname, gridI_vname,
        exgrid_fname, exgrid_vname,
        interp_style,
        elevI, maskI):

        elevI = elevI.reshape(-1)
        maskI = maskI.reshape(-1)
        cicebin.GCMRegridder_add_sheet(&self.cself,
            name.encode(),
            gridI_fname.encode(), gridI_vname.encode(),
            exgrid_fname.encode(), exgrid_vname.encode(),
            interp_style.encode(),
            <PyObject *>elevI, <PyObject *>maskI)   # Borrowed references

    def regrid_matrices(self, str sheet_name):
        cdef cicebin.RegridMatrices *crm = new cicebin.RegridMatrices(
            self.cself.ice_regridder(sheet_name.encode()))
        rm = RegridMatrices()
        rm.cself = crm
        return rm

def coo_multiply(M, xx, double fill=np.nan, ignore_nan=False):
    xx = xx.reshape(-1)
    yy = np.zeros(M._shape[0])
    yy[:] = fill

    cicebin.coo_matvec(<PyObject *>yy, <PyObject *>xx, ignore_nan,
        M._shape[0], M._shape[1],
        <PyObject *>M.row, <PyObject *>M.col, <PyObject *>M.data)

    return yy
