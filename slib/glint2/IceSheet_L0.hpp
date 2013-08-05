#pragma once

#include "IceSheet.hpp"

namespace glint2 {

class IceSheet_L0 : public IceSheet
{		// For ice model with level-value grid cells
public:
	/** Number of grid cells in the ice grid */
	size_t n2() const { return exgrid->grid2_ncells_full; }

	/** Number of grid cells in the exchange grid */
	size_t n4() const { return exgrid->ncells_full(); }

protected:

	/** Tells whether a cell in the exchange grid is masked out or not */
	bool masked(giss::HashDict<int, Cell>::iterator const &it);
	bool masked(giss::HashDict<int, Cell>::const_iterator const &it);

	/** Builds an interpolation matrix to go from height points to ice/exchange grid.
	@param overlap_type Controls matrix output to ice or exchange grid. */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_ice(IceExch dest = IceExch::ICE);

public:

	/** Adds up the (ice-covered) area of each GCM grid cell */
	virtual void accum_areas(
		giss::SparseAccumulator<int,double> &area1_m);

	/** Converts vector from ice grid (n2) to exchange grid (n4).
	The transformation is easy, and no matrix needs to be computed.
	NOTE: This only makes sense for L0 grids. */
	virtual blitz::Array<double,1> ice_to_exch(blitz::Array<double,1> const &f2);

	// exch_to_ice is not hard, but has not been needed yet.
//	virtual std::unique_ptr<giss::VectorSparseMatrix> exch_to_ice();

	/** Computes matrix to go from height-point space [nhp * n1] to atmosphere grid [n1]
	@param area1_m IN/OUT: Area of each GCM cell covered by
		(non-masked-out) ice sheet. */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_projatm(
		giss::SparseAccumulator<int,double> &area1_m);

	virtual std::unique_ptr<giss::VectorSparseMatrix> ice_to_projatm(
		giss::SparseAccumulator<int,double> &area1_m,
		IceExch src = IceExch::ICE);


	virtual boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;

};


}	// namespace glint2
