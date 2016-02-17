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

inline PyObject *RegridMatrices_regrid(RegridMatrices *cself, std::string const &spec_name, bool scale)
{
	std::unique_ptr<WeightedSparse> Mw(cself->regrid(spec_name, scale));

	PyObject *weight_py = ibmisc::cython::blitz_to_np(Mw->weight.to_dense(0));
	PyObject *M_py = ibmisc::cython::spsparse_to_tuple(Mw->M);
	return Py_BuildValue("OO", M_py, weight_py);
}

void coo_matvec(PyObject *yy_py, PyObject *xx_py, bool ignore_nan,
	size_t M_nrow, size_t M_ncol, PyObject *M_row_py, PyObject *M_col_py, PyObject *M_data_py);


}}

