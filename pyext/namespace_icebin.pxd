from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
cimport namespace_std as cstd
cimport namespace_ibmisc as cibmisc

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

cdef extern from "icebin/GCMRegridder.hpp" namespace "icebin":
	pass

	cdef cppclass GCMRegridder:
		GCMRegridder() except +
		void ncio(cibmisc.NcIO &, string vname)

cdef extern from "icebin_cython.hpp" namespace "icebin":
	cdef void GCMRegridder_init(
		GCMRegridder *self,
		string &gridA_fname,
		string &gridA_vname,
		vector[double] &hpdefs,
		bool correctA)

