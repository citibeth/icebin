#pragma once

#include <memory>
#include <glint2/Grid.hpp>
#include <glint2/matrix_ops.hpp>
#include <giss/Proj2.hpp>

/*
Height Points (HP) [nhc*n1] --A--> Ice [n2] --B--> Height Classes (HC) [nhc*n1] --C--> GCM

A: See hp2ice() below
B: See ice2hc() below
C: See get_fhc() below, then multiply by FHC in ModelE code.

*/

namespace glint2 {

/** Generates the matrices required in the GCM */
class MatrixMaker {

public:
	/** These are all left public because someone will probably want
	to look at / use them. */

	// ------------ Stuff we're passed in
	std::shared_ptr<glint2::Grid> grid1;		/// GCM Grid
	std::shared_ptr<glint2::Grid> grid2;		/// Ice Grid
//	giss::Proj2 proj;							/// GCM -> Ice Projection
	std::shared_ptr<glint2::Grid> exgrid;	/// Exchange grid (between GCM and Ice)

	std::shared_ptr<blitz::Array<int,1>> mask1;
	std::shared_ptr<blitz::Array<int,1>> mask2;

	/** Elevation of each cell Ma(L0) or vertex (L1) in the ice model */
	blitz::Array<double,1> elev2;	// [n2]

	/** Upper boundary of each height class (a, b] */
	blitz::Array<double,1> hcmax;	// [nhc-1]

	/** @return The height class of the given point in gridcell #index1. */
	int get_hc(int index1, double elev);

	/** Position of height points in elevation space (same for all GCM grid cells) */
	std::vector<double> hpdefs;	// [nhc]

	int nhc() { return hpdefs.size(); }

	/** The overlap matrix, derived from exgrid.
	NOTE: This sparse matrix has its elements in the same order as in exgrid. */
	std::unique_ptr<giss::VectorSparseMatrix> overlap_raw;	/// Overlap matrix between grid1 and grid2

	/** Masked with mask1 and mask2. */
	std::unique_ptr<giss::VectorSparseMatrix> overlap_m;

	// ------------------------------------------------

	virtual void realize();
	virtual ~MatrixMaker() {}

	// ------------------------------------------------

	virtual void compute_fhc(
		blitz::Array<double,2> *fhc1h,
		blitz::Array<double,2> *elev1h,
		blitz::Array<double,1> *fgice1) = 0;	// Portion of gridcell covered in ground ice (from landmask)

	/** Make matrix to go from
		height points [nhc*n1] to ice grid [n2].
	This is used on each ice timestep to generate SMB. */
	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_ice() = 0;

	/** Make matrix to go from
		ide grid [n2] to height classes [nhc*n1].
	The matrix hp2ice * ice2hc is used every GCM timestep to generate
	SMB for the atmosphere.  (It is later multiplied by FHC). */
	virtual std::unique_ptr<giss::VectorSparseMatrix> ice_to_hc() = 0;

	virtual std::unique_ptr<giss::VectorSparseMatrix> hp_to_hc();

	boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;

};	// class MatrixMaker

}	// namespace glint2

