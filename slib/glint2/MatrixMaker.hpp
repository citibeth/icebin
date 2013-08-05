#pragma once

#include <memory>
#include <giss/CooVector.hpp>
#include <glint2/Grid.hpp>
#include <glint2/matrix_ops.hpp>
#include <glint2/IceSheet.hpp>
#include <glint2/GridDomain.hpp>
#include <giss/hash.hpp>
#include <glint2/HCIndex.hpp>

/*
Height Points (HP) [nhc*n1] --A--> Ice [n2] --B--> Height Classes (HC) [nhc*n1] --C--> GCM

A: See hp2ice() below
B: See ice2hc() below
C: See get_fhc() below, then multiply by FHC in ModelE code.

*/

namespace glint2 {

typedef giss::SparseAccumulator<std::pair<int,int>, double, giss::HashPair<int,int>> SparseAccumulator1hc;

BOOST_ENUM_VALUES( QPAlgorithm, int,
	(SINGLE_QP)		(0)
	(MULTI_QP)		(1)
)

/** Generates the matrices required in the GCM */
class MatrixMaker
{
public:
//	std::map<int, std::unique_ptr<IceSheet>> sheets;

//	std::vector<string> sheet_names;	// Gives numbering of ALL sheets as well as names
    giss::MapDict<std::string, IceSheet> sheets;
	std::map<int, IceSheet *> sheets_by_id;
	std::unique_ptr<HCIndex> hc_index;	// Methods to extract i1 and ihc from an elevation index
protected:
	int _next_sheet_index;
	std::unique_ptr<GridDomain> domain;
public:
	bool const correct_area1;		/// Should we correct for projection and geometric error?
	HCIndex::Type _hptype;
	MatrixMaker(
		bool _correct_area1,
		std::unique_ptr<GridDomain> &&_domain)
		: _next_sheet_index(0),
		_hptype(HCIndex::Type::UNKNOWN),
		domain(std::move(_domain)),
		correct_area1(_correct_area1) {}

//	std::vector<std::unique_ptr<IceSheet>> sheets;
	
	/** These are all left public because someone will probably want
	to look at / use them. */

	/** @return Vector of names of the ice sheets. */
	std::vector<std::string> get_sheet_names() {
		std::vector<std::string> ret;
		for (auto ii=sheets.begin(); ii != sheets.end(); ++ii)
			ret.push_back(ii.key());
		return ret;
	}

	// ------------ Stuff we're passed in
	std::shared_ptr<glint2::Grid> grid1;		/// GCM Grid

	// Blitz++ arrays are reference counted internally
	// TODO: This should cover only GCM grid cells that live on our local MPI node.
	std::unique_ptr<blitz::Array<int,1>> mask1;

	/** Position of height points in elevation space (same for all GCM grid cells) */
	std::vector<double> hpdefs;	// [nhp]

	// ------------------------------------------------------
	void clear();
	int add_ice_sheet(std::unique_ptr<IceSheet> &&sheet);

	/** Call this after you've set everything up, added all ice sheets, etc. */
	void realize();

	IceSheet *operator[](int const ix) {
		auto ii = sheets_by_id.find(ix);
		if (ii == sheets_by_id.end()) return NULL;
		return ii->second;
	}
	IceSheet *operator[](std::string const &name)
		{ return sheets[name]; }

//	int nhp() const { return hpdefs.size(); }
	int n1() const { return grid1->ndata(); }
	int n3() const { return n1() * hpdefs.size(); }

	/** @return Number of elevation points for a given grid cell */
	int nhp(int i1) const { return hpdefs.size(); }
//	std::vector<int> get_used1();

	/** NOTE: Allows for two ice sheets overlapping the same GCM grid cell.
	Ice sheets cannot overlap each other (although their grids can, if we're
	guaranteed that ice-filled grid cells will never overlap). */
	void fgice(giss::CooVector<int,double> &fgice1);

	std::unique_ptr<giss::VectorSparseMatrix> hp_to_atm();

	/** @params f2 Some field on each ice grid (referenced by ID)
	TODO: This only works on one ice sheet.  Will need to be extended
	for multiple ice sheets. */
	giss::CooVector<int, double> ice_to_hp(
		std::map<int, blitz::Array<double,1>> &f2s,
		blitz::Array<double,1> &initial3,
		IceExch src = IceExch::ICE,
		QPAlgorithm qp_algorithm = QPAlgorithm::SINGLE_QP);



	boost::function<void ()> netcdf_define(NcFile &nc, std::string const &vname) const;

	virtual void read_from_netcdf(NcFile &nc, std::string const &vname);
};

std::unique_ptr<IceSheet> new_ice_sheet(Grid::Parameterization parameterization);


}	// namespace glint2

