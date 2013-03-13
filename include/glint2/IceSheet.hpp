#pragma once

#include <memory>
#include <glint2/Grid.hpp>
#include <blitz/array.h>
#include <giss/SparseMatrix.hpp>
#include <glint2/ExchangeGrid.hpp>
#include <giss/SparseAccumulator.hpp>

namespace glint2 {

class MatrixMaker;
class IceCoupler;

class IceSheet {
protected:
	friend class MatrixMaker;

	MatrixMaker *gcm;

public:
	int index;

	std::shared_ptr<glint2::Grid> grid2;		/// Ice Grid
	std::shared_ptr<glint2::ExchangeGrid> exgrid;	/// Exchange grid (between GCM and Ice)

	std::string name; //const &name() const { return grid2->name; }

	/** TODO: How does mask2 work for L1 grids? */
	std::unique_ptr<blitz::Array<int,1>> mask2;

	/** Elevation of each cell (L0) or vertex (L1) in the ice model */
	blitz::Array<double,1> elev2;	// [n2]

	// ============== Volatile / Compute Variables

	/** The overlap matrix, derived from exgrid.
	NOTE: This sparse matrix has its elements in the same order as in exgrid. */
	std::unique_ptr<giss::VectorSparseMatrix> overlap_raw;	/// Overlap matrix between grid1 and grid2

	/** Masked with mask1 and mask2. */
	std::unique_ptr<giss::VectorSparseMatrix> overlap_m;

	/** Multiply by this each ice timestep to compute SMB */
	std::unique_ptr<giss::VectorSparseMatrix> _hp_to_ice;

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

	virtual void accum_areas(
		giss::SparseAccumulator<int,double> &area1_m,
		giss::SparseAccumulator<int,double> &area1_m_hc) = 0;

	/** Make matrix to go from
		height points [nhc*n1] to ice grid [n2].
	This is used on each ice timestep to generate SMB. */
	virtual std::unique_ptr<giss::VectorSparseMatrix> compute_hp_to_ice() = 0;

	giss::VectorSparseMatrix &hp_to_ice() {
		if (!_hp_to_ice.get()) _hp_to_ice = compute_hp_to_ice();
		return *_hp_to_ice;
	}

protected:
	/** Make matrix to go from
		ide grid [n2] to height classes [nhc*n1].
	The matrix hp2ice * ice2hc is used every GCM timestep to generate
	SMB for the atmosphere.  (It is later multiplied by FHC). */
	virtual std::unique_ptr<giss::VectorSparseMatrix> ice_to_hc(
		giss::SparseAccumulator<int,double> &area1_m_hc) = 0;
public:

	virtual boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;
	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);

};	// class IceSheet

}
