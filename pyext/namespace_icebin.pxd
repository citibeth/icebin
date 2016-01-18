cimport namespace_std as cstd

# This will make the C++ class def for Rectangle available..
cdef extern from "icebin/Grid.hpp" namespace "icebin":
	pass

	cdef cppclass Vertex:
		Vertex(double, double, long) except +
		long index
		const double x
		const double y


	cdef cppclass Cell:
		Cell(cstd.vector[Vertex *]) except +

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

		
