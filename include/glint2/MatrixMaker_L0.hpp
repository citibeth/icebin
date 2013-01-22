#pragma once

#include "MatrixMaker.hpp"

namespace glint2 {

class MatrixMaker_L0 : public MatrixMaker
{		// For ice model with level-value grid cells
public:
	virtual void realize();

	/** Masked and then height classified */
	std::unique_ptr<giss::VectorSparseMatrix> overlap_m_hc;

	/** Uses: elev2, hpdefs, overlap */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_ice();

	/** Uses: elev2, hcmax, overlap */
	virtual std::unique_ptr<giss::VectorSparseMatrix> ice_to_hc();

	virtual void compute_fhc(
		blitz::Array<double,2> *fhc1h,
		blitz::Array<double,1> *fgice1);	// Portion of gridcell covered in ground ice (from landmask)

};


}	// namespace glint2
