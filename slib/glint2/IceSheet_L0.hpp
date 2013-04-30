#pragma once

#include "IceSheet.hpp"

namespace glint2 {

class IceSheet_L0 : public IceSheet
{		// For ice model with level-value grid cells
public:
	virtual void realize();

	long n1() const { return exgrid->grid1_ncells_full; }
	long n2() const { return exgrid->grid2_ncells_full; }
	long n3() const { return exgrid->ncells_full(); }

protected:

	/** Tells whether a cell in the exchange grid is masked out or not */
	bool masked(HashDict<int, Cell>::iterator const &it)
	bool masked(HashDict<int, Cell>::const_iterator const &it)

	/** Builds an interpolation matrix to go from height points to ice/exchange grid.
	@param overlap_type Controls matrix output to ice or exchange grid. */
	std::unique_ptr<giss::VectorSparseMatrix> hp_interp(Overlap overlap_type)

public:

	/** Adds up the (ice-covered) area of each GCM grid cell */
	virtual void accum_areas(
		giss::SparseAccumulator<int,double> &area1_m);

	/** Computes matrix to go from height-point space [nhp * n1] to ice grid [n2] */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_ice();

	/** Computes matrix to go from height-point space [nhp * n1] to atmosphere grid [n1]
	@param area1_m IN/OUT: Area of each GCM cell covered by
		(non-masked-out) ice sheet. */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_atm(
		giss::SparseAccumulator<int,double> &area1_m);




	virtual boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;

};


}	// namespace glint2
