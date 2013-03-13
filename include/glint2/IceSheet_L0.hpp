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


	/**
	@param area1_m IN/OUT: Area of each GCM cell covered by
		(non-masked-out) ice sheet.
	@param area1_m_hc IN/OUT: Area of each GCM cell / height class coverd by
		(non-masked-out) ice sheet
	    NOTE: Indexed in 1-D according to HCIndex convention [nhc * n1] */
	void accum_areas(
		giss::SparseAccumulator<int,double> &area1_m,
		giss::SparseAccumulator<int,double> &area1_m_hc);

	virtual boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;

};


}	// namespace glint2
