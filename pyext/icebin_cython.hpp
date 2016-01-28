#pragma once

#include <icebin/GCMRegridder.hpp>

namespace icebin {

inline void GCMRegridder_init(GCMRegridder *self,
	std::string const &gridA_fname,
	std::string const &gridA_vname,
	std::vector<double> &hpdefs,
	bool _correctA)
{
	// Read gridA
	ibmisc::NcIO ncio(gridA_fname, netCDF::NcFile::read);
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
		ibmisc::Domain<int>(
			std::vector<int>(gridA->indexing.base), std::move(high)),
		std::move(hpdefs),
		ibmisc::Indexing<long,long>(
			{0,0}, {gridA->ndata(), nhp}, {1,0}),
		_correctA);

	self->hpdefs = hpdefs;
}

};
