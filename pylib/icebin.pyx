# Get definitions of C++ classes in the cpp namespace.
#cimport cpp
cimport namespace_std as cstd
cimport namespace_icebin as cicebin
from cython.operator cimport dereference as deref, preincrement as inc
cimport ibmisc

cdef class GCMRegridder:
	cdef cicebin.GCMRegridder cself

	def __init__(self, gridA_fname, gridA_vname, hpdefs, correctA):
		cicebin.GCMRegridder_init(&self.cself, gridA_fname.encode(), gridA_vname.encode(), hpdefs, correctA)

	def ncio(self, ibmisc.NcIO ncio, vname):
		self.cself.ncio(deref(ncio.cself), vname.encode())



# cpdef test_icebin():
# 	cdef cicebin.Grid grid
# 	cdef cicebin.Vertex *vertex
# 	cstd.vector[cicebin.Vertex *] vertices	# For one cell
# 
# 	vcache = dict()
# 	for i in range(0,2):
# 		x = i*2.
# 		for j in range(0,2):
# 			j*1.
# 
# 			vertices.clear()
# 			
# 
# 			vertex = vcache[(i,j)] if (i,j) in vcache else 
# 			if (i,j) in vcache:
# 
# 
# 	cdef cicebin.Vertex *vertex;
# 	vertex = grid.vertices.add_claim(new cicebin.Vertex(1.,2.,-1))
# 	print(vertex.x)
# 
# #	cdef cicebin.Vertex vertex(1.,2.,-1)
# #	vertex.x = 1.
# #	vertex.y = 2.
# #	grid.vertices.add(vertex)
# #	grid.vertices.add(cicebin.Vertex(1.,2.,-1))
# #	cdef cicebin.Vertex *vertex = grid.vertices.at(0)
# 