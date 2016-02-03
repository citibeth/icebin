cimport cicebin
cimport cibmisc		# C++ stuff in ibmisc
cimport ibmisc		# Cython stuff in ibmisc
import numpy as np
import scipy.sparse

from cython.operator cimport dereference as deref, preincrement as inc

cdef class IceRegridder:
	pass

cdef class RegridMatrices:
	cdef cicebin.RegridMatrices *cself

	def __dealloc__(self):
		del self.cself

	def regrid(self, spec_name):
		data,shape = cicebin.RegridMatrices_regrid(self.cself, spec_name.encode())
		# scipy.sparse.coo_matrix((data1, (rows1, cols1)), shape=(nrow1, ncol1))
		return scipy.sparse.coo_matrix(data, shape)

	def weight(self, spec_name, fill_value=0):
		data,shape = cicebin.RegridMatrices_weight(self.cself, spec_name.encode(), fill_value)
		ret = np.zeros((shape[0],))
		ret[:] = fill_value
		cdef int ix
		cdef double val
		for ix,val in zip(data[1][0], data[0]):
			ret[ix] = val
		return ret

cdef class GCMRegridder:
	cdef cicebin.GCMRegridder cself

	def __init__(self, gridA_fname, gridA_vname, hpdefs, correctA):
		cicebin.GCMRegridder_init(&self.cself, gridA_fname.encode(), gridA_vname.encode(), hpdefs, correctA)

	def ncio(self, ibmisc.NcIO ncio, vname):
		self.cself.ncio(deref(ncio.cself), vname.encode())

	def add_sheet(self, name,
		gridI_fname, gridI_vname,
		exgrid_fname, exgrid_vname,
		interp_style,
		elevI):
		pass

	def regrid_matrices(self, str sheet_name):
		cdef cicebin.RegridMatrices *crm = new cicebin.RegridMatrices(
			self.cself.sheet(sheet_name.encode()))
		rm = RegridMatrices()
		rm.cself = crm
		return rm

