# On reference counting...
# https://groups.google.com/forum/#!topic/cython-users/pZFurj3rUyg

from cpython.object cimport *		# PyObject
cimport cicebin
cimport cibmisc		# C++ stuff in ibmisc
cimport ibmisc		# Cython stuff in ibmisc
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

	def regrid(self, str spec_name):
		data,shape = cicebin.RegridMatrices_regrid(self.cself, spec_name.encode())
		# scipy.sparse.coo_matrix((data1, (rows1, cols1)), shape=(nrow1, ncol1))
		return scipy.sparse.coo_matrix(data, shape)

	def weight(self, spec_name, fill_value=0):
		data,shape = cicebin.RegridMatrices_scale(self.cself, spec_name.encode(), fill_value)
		ret = np.zeros((shape[0],))
		ret[:] = fill_value
		cdef int ix
		cdef double val
		for ix,val in zip(data[1][0], data[0]):
			ret[ix] = val
		return ret

cdef class GCMRegridder:
	cdef cicebin.GCMRegridder cself

	def __init__(self, *args):
		cdef ibmisc.NcIO ncio

		try:
			gridA_fname, gridA_vname, hpdefs, correctA = args
			cicebin.GCMRegridder_init(&self.cself, gridA_fname.encode(), gridA_vname.encode(), hpdefs, correctA)
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
			<PyObject *>elevI, <PyObject *>maskI)	# Borrowed references

	def regrid_matrices(self, str sheet_name):
		cdef cicebin.RegridMatrices *crm = new cicebin.RegridMatrices(
			self.cself.sheet(sheet_name.encode()))
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
