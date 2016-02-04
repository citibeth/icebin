from cpython.object cimport *	# PyObject
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
cimport cibmisc

# This will make the C++ class def for Rectangle available..
cdef extern from "icebin/Grid.hpp" namespace "icebin":
	pass

	cdef cppclass Vertex:
		Vertex(double, double, long) except +
		long index
		const double x
		const double y


	cdef cppclass Cell:
		Cell(vector[Vertex *]) except +

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

	cdef cppclass IceRegridder:
		pass

	cdef cppclass GCMRegridder:
		GCMRegridder() except +
		void ncio(cibmisc.NcIO &, string vname) except +
		IceRegridder *sheet(string name) except +


	cdef cppclass RegridMatrices:
		RegridMatrices(IceRegridder *) except +


cdef extern from "icebin_cython.hpp" namespace "icebin::cython":
	cdef void GCMRegridder_init(
		GCMRegridder *self,
		string &gridA_fname,
		string &gridA_vname,
		vector[double] &hpdefs,
		bool correctA) except +

	cdef void GCMRegridder_add_sheet(
		GCMRegridder *cself,
		string &name,
		string &gridI_fname, string &gridI_vname,
		string &exgrid_fname, string &exgrid_vname,
		string &sinterp_style,
		PyObject *elevI_py, PyObject *maskI_py) except +		# PyObject=Borrowed reference, object = owned reference

	cdef object RegridMatrices_regrid(RegridMatrices *self, string spec_name) except +

	cdef object RegridMatrices_weight(RegridMatrices *self, string spec_name, double fill_value) except +
