cimport cicebin
cimport cibmisc		# C++ stuff in ibmisc
cimport ibmisc		# Cython stuff in ibmisc

from cython.operator cimport dereference as deref, preincrement as inc

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

