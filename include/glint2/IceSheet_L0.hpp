#pragma once

#include "IceSheet.hpp"

namespace glint2 {

class IceSheet_L0 : public IceSheet
{		// For ice model with level-value grid cells
public:
	virtual void realize();

	long n1() const { return exgrid->grid1_ncells_full; }
	long n2() const { return exgrid->grid2_ncells_full; }

	/** Masked and then height classified */
	std::unique_ptr<giss::VectorSparseMatrix> overlap_m_hc;

	/** Uses: elev2, hpdefs, overlap */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_ice();

	/** Uses: elev2, hcmax, overlap */
	virtual std::unique_ptr<giss::VectorSparseMatrix> ice_to_hc();

	virtual void compute_fhc(
		blitz::Array<double,2> *fhc1h,
		blitz::Array<double,1> *fgice1);	// Portion of gridcell covered in ground ice (from landmask)

	virtual boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;

};


}	// namespace glint2
