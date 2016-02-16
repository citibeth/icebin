#pragma once

#include <Python.h>
#include <ibmisc/cython.hpp>
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
	PyObject *elevI_py, PyObject *maskI_py);

inline PyObject *RegridMatrices_regrid(RegridMatrices *cself, std::string const &spec_name)
{
	SparseMatrix M(cself->regrid(spec_name));
	return ibmisc::cython::spsparse_to_tuple(M);
}

inline PyObject *RegridMatrices_scale(RegridMatrices *cself, std::string const &spec_name, double fill_value)
{
	SparseVector w(cself->scale(spec_name));
	// TODO: Copy only once instead of twice.  Don't bother for now...
	return ibmisc::cython::blitz_to_np(w.to_dense(fill_value));
}

void coo_matvec(PyObject *yy_py, PyObject *xx_py, bool ignore_nan,
	size_t M_nrow, size_t M_ncol, PyObject *M_row_py, PyObject *M_col_py, PyObject *M_data_py);


}}
