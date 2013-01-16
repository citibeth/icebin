#pragma once

#include "MatrixMaker.hpp"

namespace glint2 {

class MatrixMaker_L0 : public MatrixMaker
{		// For ice model with level-value grid cells
public:
	/** Masked and then height classified */
	std::unique_ptr<VectorSparseMatrix> overlap_m_hc;

	virtual MatrixMaker_L0(MatrixMakerData &&data);

	/** Uses: elev2, hpdefs, overlap */
	virtual std::unique_ptr<VectorSparseMatrix> hp_to_ice();

	/** Uses: elev2, hcmax, overlap */
	virtual std::unique_ptr<VectorSparseMatrix> ice_to_hc();

};


}	// namespace glint2
