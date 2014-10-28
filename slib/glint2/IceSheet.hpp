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

#include <unordered_set>
#include <memory>
#include <glint2/Grid.hpp>
#include <blitz/array.h>
#include <giss/SparseMatrix.hpp>
#include <glint2/ExchangeGrid.hpp>
#include <giss/SparseAccumulator.hpp>
#include <giss/CooVector.hpp>

namespace glint2 {

enum class ProjCorrect {NATIVE_TO_PROJ=0, PROJ_TO_NATIVE=1};
enum class FactorUse {MULTIPLY=0, DIVIDE=1};

/** For subroutines that can do things to/from either the ice grid or
the interpolation grid */
BOOST_ENUM_VALUES( IceInterp, int,
	(ICE)		(0)
	(INTERP)	(1)
)

/** Controls how we interpolate from elevation class space to the ice grid */
BOOST_ENUM_VALUES( InterpStyle, int,
	(Z_INTERP)			(0)
	(ELEV_CLASS_INTERP)	(1)
	(BILIN_INTERP)		(2)
)

class MatrixMaker;
class IceCoupler;

/** NOTE: The Interpolation Grid (grid4) is the grid to which we
interpolate, and which ultimately defines the basis functions in
elevation point space.  By default, this is the Ice Grid.  But for
L0 ice grids, we can use the Exchange Grid as the Interpolation Grid,
in order to keep the RM matrix local. */
class IceSheet {
protected:
	friend class MatrixMaker;

	MatrixMaker *gcm;

public:
	int index;

	InterpStyle interp_style;

	std::shared_ptr<glint2::Grid> grid2;		/// Ice Grid
	std::shared_ptr<glint2::ExchangeGrid> exgrid;	/// Exchange grid (between GCM and Ice)

	std::string name; //const &name() const { return grid2->name; }

	/** TODO: How does mask2 work for L1 grids? */
	std::unique_ptr<blitz::Array<int,1>> mask2;

	/** Elevation of each cell (L0) or vertex (L1) in the ice model */
	blitz::Array<double,1> elev2;	// [n2]

	// ===================================================

	// -------------------------------------------

	IceSheet();
	virtual void clear();
protected:
	virtual void realize();
public:
	void filter_cells1(boost::function<bool (int)> const &include_cell1);

	virtual ~IceSheet();

	// ------------------------------------------------
	/** Number of dimensions of atmosphere vector space */
	virtual size_t n1() const;

	/** Number of dimensions of ice vector space */
	virtual size_t n2() const = 0;

	/** Number of dimensions of interpolation grid vector space. */
	virtual size_t n4() const
		{ return n2(); }

	// ------------------------------------------------
	/** Diagonal matrix converts values from native atmosphere grid to projected atmosphere grid (or vice versa)
	@param direction Direction to convert vectors (NATIVE_TO_PROJ or PROJ_TO_NATIVE) */
	std::unique_ptr<giss::VectorSparseMatrix> atm_proj_correct(ProjCorrect direction);

	/** Puts GCM grid correction factors into the area1_m variable often
	involved in regridding matrices.  Avoides having to create a new sparse
	matrix just for this purpose.  This subroutine is only really useful when
	working with one ice sheet at a time. */
	void atm_proj_correct(
		giss::SparseAccumulator<int,double> &area1_m,
		ProjCorrect direction);

	// ------------------------------------------------

	/** Adds up the (ice-covered) area of each GCM grid cell */
	virtual void accum_areas(
		giss::SparseAccumulator<int,double> &area1_m) = 0;

	/** Computes matrix to go from elevation point space [nhp * n1] to ice grid [n2] */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_iceinterp(IceInterp dest) = 0;

	/** Convert from ice grid to the interpolation grid.
	NOTE: This simple/stub implementation will generally be overridden. */
	virtual blitz::Array<double,1> const ice_to_interp(blitz::Array<double,1> const &f2)
		{ return f2; }

	/** Computes matrix to go from elevation point space [nhp * n1]
	to projected atmosphere grid [n1].  NOTE: Corrections for geometric
	and projection error when going between Cartesian and Spherical
	space are not accounted for here.
	@param area1_m IN/OUT: Area of each GCM cell covered by
		(non-masked-out) ice sheet.  Must divide result by this number. */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_projatm(
		giss::SparseAccumulator<int,double> &area1_m) = 0;

	virtual std::unique_ptr<giss::VectorSparseMatrix> iceinterp_to_projatm(
		giss::SparseAccumulator<int,double> &area1_m,
		IceInterp src) = 0;

public:

	virtual boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;
	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);

};	// class IceSheet

// ----------------------------------------------------

extern void linterp_1d(
	std::vector<double> const &xpoints,
	double xx,
	int *indices, double *weights);	// Size-2 arrays


/** We only really expect this to work for Greenland.  Don't worry
about south pole in lon/lat coordinates and Antarctica.
[n2 x (nhc * n1)] sparse matrix */
extern std::unique_ptr<giss::VectorSparseMatrix> 
bilin_interp(
MatrixMaker *gcm,
Grid const &grid1_lonlat,
Grid const &grid2,
std::vector<double> const &hpdefs,
blitz::Array<double,1> const &elev2,
blitz::Array<int,1> const *mask1,		// [n1] Shows where we will / will not expect landice
blitz::Array<int,1> const *mask2);


}
