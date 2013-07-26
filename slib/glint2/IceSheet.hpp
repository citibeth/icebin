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

enum class ProjCorrect {NATIVE_TO_PROJ, PROJ_TO_NATIVE};

class MatrixMaker;
class IceCoupler;

class IceSheet {
protected:
	friend class MatrixMaker;

	MatrixMaker *gcm;

	enum class Overlap {ICE, EXCH};

public:
	int index;

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

	// ------------------------------------------------
	/** Diagonal matrix converts values from native atmosphere grid to projected atmosphere grid (or vice versa)
	@param direction Direction to convert vectors (NATIVE_TO_PROJ or PROJ_TO_NATIVE) */
	std::unique_ptr<giss::VectorSparseMatrix> atm_proj_correct(ProjCorrect direction);

	// ------------------------------------------------

	/** Adds up the (ice-covered) area of each GCM grid cell */
	virtual void accum_areas(
		giss::SparseAccumulator<int,double> &area1_m) = 0;

	/** Computes matrix to go from height-point space [nhp * n1] to ice grid [n2] */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_ice() = 0;

	/** Computes matrix to go from height-point space [nhp * n1]
	to projected atmosphere grid [n1].  NOTE: Corrections for geometric
	and projection error when going between Cartesian and Spherical
	space are not accounted for here.
	@param area1_m IN/OUT: Area of each GCM cell covered by
		(non-masked-out) ice sheet.  Must divide result by this number. */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_projatm(
		giss::SparseAccumulator<int,double> &area1_m) = 0;

	virtual std::unique_ptr<giss::VectorSparseMatrix> ice_to_projatm(
		giss::SparseAccumulator<int,double> &area1_m) = 0;

public:

	virtual boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;
	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);

};	// class IceSheet

}
