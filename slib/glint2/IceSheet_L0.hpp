/*
 * GLINT2: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013 by Robert Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
	@param dest Controls matrix output to ice or exchange grid.
	@param fill_masked If true, then ice/exch grid cells that are masked out will be treated as
		if they have an elevation point of -1.  This can be used later, where the matrix is applied,
		to fill in a background field for masked ice grid cells.
	*/
	std::unique_ptr<giss::VectorSparseMatrix> hp_to_iceexch(IceExch dest, bool fill_masked = false);

public :
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_iceinterp(
		IceInterp dest, bool fill_masked)
	{
		IceExch iedest = (dest == IceInterp::ICE ? IceExch::ICE : interp_grid);
		return hp_to_iceexch(iedest, fill_masked);
	}


public:

	/** Adds up the (ice-covered) area of each GCM grid cell */
	virtual void accum_areas(
		giss::MapSparseVector<int,double> &area1_m);

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
		giss::MapSparseVector<int,double> &area1_m);

protected :
	std::unique_ptr<giss::VectorSparseMatrix> iceexch_to_projatm(
		giss::MapSparseVector<int,double> &area1_m,
		IceExch src = IceExch::ICE);

public:
	virtual std::unique_ptr<giss::VectorSparseMatrix> iceinterp_to_projatm(
		giss::MapSparseVector<int,double> &area1_m,
		IceInterp src)
	{
		IceExch iesrc = (src == IceInterp::ICE ? IceExch::ICE : interp_grid);
		return iceexch_to_projatm(area1_m, iesrc);
	}


	virtual boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;
	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);

};


}	// namespace glint2
