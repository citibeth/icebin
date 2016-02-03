#pragma once

#include <ibmisc/netcdf.hpp>
#include <ibmisc/cython.hpp>
#include "icebin_cython.hpp"

using namespace ibmisc;
using namespace ibmisc::cython;

namespace icebin {
namespace cython {

void GCMRegridder_init(GCMRegridder *self,
	std::string const &gridA_fname,
	std::string const &gridA_vname,
	std::vector<double> &hpdefs,
	bool _correctA)
{
	// Read gridA
	NcIO ncio(gridA_fname, netCDF::NcFile::read);
	std::unique_ptr<Grid> gridA = read_grid(ncio, gridA_vname);
	ncio.close();

	// Construct a domain to be the full extent of indexing.
	int rank = gridA->indexing.base.size();
	std::vector<int> high(rank);
	for (int k=0; k<rank; ++k)
		high[k] = gridA->indexing.base[k] + gridA->indexing.extent[k];


	// Put it together
	long nhp = hpdefs.size();
	self->init(
		std::move(gridA),
		Domain<int>(
			std::vector<int>(gridA->indexing.base), std::move(high)),
		std::move(hpdefs),
		Indexing<long,long>(
			{0,0}, {gridA->ndata(), nhp}, {1,0}),
		_correctA);

	self->hpdefs = hpdefs;
}

void GCMRegridder_add_sheet(GCMRegridder *self,
	std::string name,
	std::string const &gridI_fname, std::string const &gridI_vname,
	std::string const &exgrid_fname, std::string const &exgrid_vname,
	std::string const &sinterp_style,
	PyObject *elevI_py)
{
	NcIO ncio_I(gridI_fname, netCDF::NcFile::read);
	std::unique_ptr<Grid> gridI(read_grid(ncio_I, "grid"));
	ncio_I.close();

	NcIO ncio_exgrid(exgrid_fname, netCDF::NcFile::read);
	std::unique_ptr<Grid> exgrid(read_grid(ncio_exgrid, "grid"));
	ncio_exgrid.close();

	auto interp_style(parse_enum<InterpStyle>(sinterp_style));
	SparseVector elevI_sp;
	spsparse::to_sparse(elevI_sp,
		np_to_blitz<double,1>(elevI_py, "elevI", {gridI->ndata()}));

	auto sheet(new_ice_regridder(gridI->parameterization));
	sheet->init(name, std::move(gridI), std::move(exgrid),
		interp_style, std::move(elevI_sp));
	self->add_sheet(std::move(sheet));
}



}}

