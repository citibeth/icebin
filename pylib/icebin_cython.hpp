#pragma once

#include <Python.h>
#include <icebin/GCMRegridder.hpp>

namespace icebin {
namespace cython {

extern void GCMRegridder_init(GCMRegridder *self,
	std::string const &gridA_fname,
	std::string const &gridA_vname,
	std::vector<double> &hpdefs,
	bool _correctA);

extern void GCMRegridder_add_sheet(GCMRegridder *self,
	std::string name,
	std::string const &gridI_fname, std::string const &gridI_vname,
	std::string const &exgrid_fname, std::string const &exgrid_vname,
	std::string const &sinterp_style,
	PyObject *elevI_py);


}}

