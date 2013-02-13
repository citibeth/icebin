#pragma once

#include <memory>
#include <glint2/Grid.hpp>
#include <glint2/matrix_ops.hpp>
#include <glint2/IceSheet.hpp>

/*
Height Points (HP) [nhc*n1] --A--> Ice [n2] --B--> Height Classes (HC) [nhc*n1] --C--> GCM

A: See hp2ice() below
B: See ice2hc() below
C: See get_fhc() below, then multiply by FHC in ModelE code.

*/

namespace glint2 {


/** Generates the matrices required in the GCM */
class MatrixMaker
{
public:
	std::vector<std::unique_ptr<IceSheet>> sheets;
	
	/** These are all left public because someone will probably want
	to look at / use them. */

	// ------------ Stuff we're passed in
	std::shared_ptr<glint2::Grid> grid1;		/// GCM Grid

	// Blitz++ arrays are reference counted internally
	std::unique_ptr<blitz::Array<int,1>> mask1;

	/** Position of height points in elevation space (same for all GCM grid cells) */
	std::vector<double> hpdefs;	// [nhc]

	/** Upper boundary of each height class (a, b] */
	blitz::Array<double,1> hcmax;	// [nhc-1]

	// ------------------------------------------------------
	void clear();
	int add_ice_sheet(std::unique_ptr<IceSheet> &&sheet);

	/** Call this after you've set everything up, added all ice sheets, etc. */
	void realize();

	IceSheet &operator[](int ix) { return *sheets[ix]; }
	IceSheet &operator[](std::string const &name) {
		for (auto ii=sheets.begin(); ii != sheets.end(); ++ii) {
			if (name == (*ii)->name) return **ii;
		}
		fprintf(stderr, "Could not find ice sheet named: %s\n", name.c_str());
		throw std::exception();
	}

	/** @return The height class of the given point in gridcell #index1. */
	int get_hc(int index1, double elev);

	int nhc() { return hpdefs.size(); }
	int n1() { return grid1->ndata(); }

	/** NOTE: Does not necessarily assume that ice sheets do not overlap on the same GCM grid cell */
	void compute_fhc(
		blitz::Array<double,2> *fhc1h,	// OUT
		blitz::Array<double,1> *fgice1);	// OUT: Portion of gridcell covered in ground ice (from landmask)


	boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;

	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);
};

std::unique_ptr<IceSheet> new_ice_sheet(Grid::Parameterization parameterization);


}	// namespace glint2

