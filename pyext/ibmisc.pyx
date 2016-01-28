cimport namespace_ibmisc as cibmisc
cimport namespace_std as cstd
from cython.operator cimport dereference as deref, preincrement as inc

cdef class NcIO:
	cdef cibmisc.NcIO *cself

	def __cinit__(self, filePath, sfMode):
		self.cself = cibmisc.new_ncio(filePath.encode(), sfMode.encode());

	def __dealloc__(self):
		del self.cself

	def close(self):
		self.cself.close()
