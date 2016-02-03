#pragma once

#include <Python.h>
#include <icebin/GCMRegridder.hpp>

namespace icebin {
namespace cython {

extern void GCMRegridder_init(GCMRegridder *cself,
	std::string const &gridA_fname,
	std::string const &gridA_vname,
	std::vector<double> &hpdefs,
	bool _correctA);

extern void GCMRegridder_add_sheet(GCMRegridder *cself,
	std::string name,
	std::string const &gridI_fname, std::string const &gridI_vname,
	std::string const &exgrid_fname, std::string const &exgrid_vname,
	std::string const &sinterp_style,
	PyObject *elevI_py);

inline Pyobject *RegridMatrices_regrid(RegridMatrices *cself, std::string const &spec_name)
{
	SparseMatrix M;
	cself->regrid(M, spec_name);
	return spsparse_to_tuple(M);
}

inline Pyobject *RegridMatrices_weight(RegridMatrices *cself, std::string const &spec_name, double fill_value)
{
	SparseVector w;
	cself->weight(w, spec_name);
	// TODO: Copy only once instead of twice.  Don't bother for now...
	return blitz_to_np(w.to_dense(fill_value))
}


}}

