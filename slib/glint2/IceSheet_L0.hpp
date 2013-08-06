#pragma once

#include "IceSheet.hpp"

namespace glint2 {

BOOST_ENUM_VALUES( IceExch, int,
	(ICE)	(0)
	(EXCH)	(1)
)

class IceSheet_L0 : public IceSheet
{		// For ice model with level-value grid cells
public:
	/** The grid we use at the interpolation grid (exchange or ice) */
	IceExch interp_grid;

	/** Number of grid cells in the ice grid */
	size_t n2() const { return exgrid->grid2_ncells_full; }

	/** Number of grid cells in the interpolation grid */
	size_t n4() const
		{ return interp_grid == IceExch::ICE ? n2() : exgrid->ncells_full(); }

	IceSheet_L0() : interp_grid(IceExch::EXCH) {}

protected:

	/** Number of grid cells in the exchange grid */
	size_t niceexch(IceExch grid) const
		{ return grid == IceExch::ICE ? n2() : exgrid->ncells_full(); }

	/** Tells whether a cell in the exchange grid is masked out or not */
	bool masked(giss::HashDict<int, Cell>::iterator const &it);
	bool masked(giss::HashDict<int, Cell>::const_iterator const &it);

protected :
	/** Builds an interpolation matrix to go from height points to ice/exchange grid.
	@param overlap_type Controls matrix output to ice or exchange grid. */
	std::unique_ptr<giss::VectorSparseMatrix> hp_to_iceexch(IceExch dest);

public :
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_iceinterp(
		IceInterp dest)
	{
		IceExch iedest = (dest == IceInterp::ICE ? IceExch::ICE : interp_grid);
printf("hp_to_iceinterp(): dest=%s, iedest=%s\n", dest.str(), iedest.str());
		return hp_to_iceexch(iedest);
	}


public:

	/** Adds up the (ice-covered) area of each GCM grid cell */
	virtual void accum_areas(
		giss::SparseAccumulator<int,double> &area1_m);

	/** Converts vector from ice grid (n2) to exchange grid (n4).
	The transformation is easy, and no matrix needs to be computed.
	NOTE: This only makes sense for L0 grids. */
	blitz::Array<double,1> const ice_to_interp(blitz::Array<double,1> const &f2);

	// exch_to_ice is not hard, but has not been needed yet.
//	virtual std::unique_ptr<giss::VectorSparseMatrix> exch_to_ice();

	/** Computes matrix to go from height-point space [nhp * n1] to atmosphere grid [n1]
	@param area1_m IN/OUT: Area of each GCM cell covered by
		(non-masked-out) ice sheet. */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_projatm(
		giss::SparseAccumulator<int,double> &area1_m);

protected :
	std::unique_ptr<giss::VectorSparseMatrix> iceexch_to_projatm(
		giss::SparseAccumulator<int,double> &area1_m,
		IceExch src = IceExch::ICE);

public:
	virtual std::unique_ptr<giss::VectorSparseMatrix> iceinterp_to_projatm(
		giss::SparseAccumulator<int,double> &area1_m,
		IceInterp src)
	{
		IceExch iesrc = (src == IceInterp::ICE ? IceExch::ICE : interp_grid);
		return iceexch_to_projatm(area1_m, iesrc);
	}


	virtual boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;
	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);

};


}	// namespace glint2
