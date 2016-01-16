
# This will make the C++ class def for Rectangle available..
cdef extern from "icebin/Grid.hpp" namespace "icebin":

	# Now, let’s add the Rectangle class to this extern from block -
    # just copy the class name from Rectangle.h and adjust for Cython
    # syntax, so now it becomes:
	cdef cppclass Vertex:
		# We now need to declare the attributes and methods for use on Cython:
		Vertex(double, double, int) except +
		long index
		const double x
		const double y

cdef class PyVertex:
	cdef Vertex *thisptr
	def __cinit__(self, double x, double y, int index):
		self.thisptr = new Vertex(x,y,index)
	def __dealloc__(self):
		del self.thisptr

	@property
	def x(self):
		return self.thisptr.x

	@property
	def y(self):
		return self.thisptr.y

	@property
	def index(self):
		return self.thisptr.index

