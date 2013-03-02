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
//	std::map<int, std::unique_ptr<IceSheet>> sheets;
	giss:MapDict<int, IceSheet> sheets;
protected:
	int _next_sheet_index;
	std::unique_ptr<GridDomain> domain;
public:
	MatrixMaker(std::unique_ptr<GridDomain> &&_domain) : domain(std::move(_domain)) {}

//	std::vector<std::unique_ptr<IceSheet>> sheets;
	
	/** These are all left public because someone will probably want
	to look at / use them. */

	// ------------ Stuff we're passed in
	std::shared_ptr<glint2::Grid> grid1;		/// GCM Grid

	// Blitz++ arrays are reference counted internally
	// TODO: This should cover only GCM grid cells in our domain.
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

	/** NOTE: Allows for two ice sheets overlapping the same GCM grid cell.
	Ice sheets cannot overlap each other (although their grids can, if we're
	guaranteed that ice-filled grid cells will never overlap). */
	void compute_fhc(
		giss::CooVector<int,double> &fhc1h,
		giss::CooVector<int,double> &fgice1);


	boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;

	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);
};

std::unique_ptr<IceSheet> new_ice_sheet(Grid::Parameterization parameterization);


}	// namespace glint2

