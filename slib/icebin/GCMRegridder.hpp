#pragma once

#include <functional>
#include <unordered_set>
#include <blitz/array.h>

#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <spsparse/eigen.hpp>

#include <icebin/Grid.hpp>
#include <icebin/sparse.hpp>

namespace icebin {

/** Controls how we interpolate from elevation class space to the ice grid */
BOOST_ENUM_VALUES( InterpStyle, int,
	(Z_INTERP)			(0)
	(ELEV_CLASS_INTERP)	(1)
)

BOOST_ENUM_VALUES( Weighting, int,
	(NONE)	(0)
	(WHOLE_CELL)	(1)
	(PARTIAL_CELL)	(2)
)


class GCMRegridder;

/** Creates regridding matrices for a single ice sheet. */
class IceRegridder {
public:
	typedef Grid::Parameterization Type;

	friend class GCMRegridder;
	friend class UrA;
	friend class UrE;
	GCMRegridder const *gcm;	/// Parent pointer
protected:

	Type type;
	std::string name;	/// "greenland", "antarctica", etc.
	std::unique_ptr<Grid> gridI;			/// Ice grid outlines
	std::unique_ptr<Grid> exgrid;		/// Exchange grid outlines (between GCM and Ice)
	InterpStyle interp_style;	/// How we interpolate I<-E

	/** Elevation of grid cells in ice grid (I).
	This also implies a mask: cells not listed in this SparseVector are masked out. */
	SparseVector elevI;
	// ---------------------------------

	// Functions used by corresponding functions in GCMRegridder
	/** Remove unnecessary GCM grid cells. */
	void filter_cellsA(std::function<bool(long)> const &keepA);

public:
	IceRegridder();
	void clear();
	void init(
		std::string const &_name,
		std::unique_ptr<Grid> &&_gridI,
		std::unique_ptr<Grid> &&_exgrid,
		InterpStyle _interp_style,
		SparseVector &&elevI);

	virtual ~IceRegridder();
	std::unordered_map<long,double> elevI_hash() const;

	// ------------------------------------------------
	/** Number of dimensions of ice vector space */
	virtual size_t nI() const = 0;

	/** Number of dimensions of interpolation grid vector space. */
	virtual size_t nG() const = 0;

	// ------------------------------------------------

	/** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
	NOTE: wAvAp == sApvA */
	void sApvA(spsparse::SparseTriplets<SparseMatrix> &w, std::function<bool(long)> const &filter_fn);

	/** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
	NOTE: wAvAp == sApvA */
	void sEpvE(spsparse::SparseTriplets<SparseMatrix> &w, std::function<bool(long)> const &filter_fn);

	virtual void GvEp(
		spsparse::SparseTriplets<SparseMatrix> &ret) const = 0;

	virtual void GvI(
		spsparse::SparseTriplets<SparseMatrix> &ret) const = 0;

	virtual void GvAp(
		spsparse::SparseTriplets<SparseMatrix> &ret) const = 0;

	virtual void ncio(ibmisc::NcIO &ncio, std::string const &vname);

};	// class IceRegridder

std::unique_ptr<IceRegridder> new_ice_regridder(IceRegridder::Type type);


// ----------------------------------------------------

extern void linterp_1d(
	std::vector<double> const &xpoints,
	double xx,
	int *indices, double *weights);	// Size-2 arrays

// ================================================
/** Generates the matrices required in the GCM */
class GCMRegridder
{
public:
	std::unique_ptr<Grid> gridA;

//	ibmisc::Domain<int> domainA;				// What's in our MPI halo?
	bool correctA;		/// Should we correct for projection and geometric error?


	/** Convert between (iA, iHP) <--> (iE) */
	ibmisc::Indexing<long,long> indexingHP;
	/** Position of height points in elevation space (same for all GCM
	grid cells) */
	std::vector<double> hpdefs;	// [nhp]

	IndexSet<std::string> sheets_index;
	typedef std::vector<std::unique_ptr<IceRegridder>> SheetsT;
	SheetsT sheets;

public:
	GCMRegridder() {}

	void clear();
	void init(
		std::unique_ptr<Grid> &&_gridA,
//		ibmisc::Domain<int> &&_domainA,		// Tells us which cells in gridA to keep...
		std::vector<double> &&_hpdefs,
		ibmisc::Indexing<long,long> &&_indexingHP,
		bool _correctA);

	// -----------------------------------------

	void add_sheet(std::unique_ptr<IceRegridder> &&sheet)
	{
		printf("Adding IceRegridder: '%s'\n", sheet->name.c_str());
		sheet->gcm = this;
		size_t ix = sheets_index.add(sheet->name);
		sheets.push_back(std::move(sheet)));
	}

	void add_sheet(std::string name, std::unique_ptr<IceRegridder> &&sheet)
	{
		sheet->name = name;
		add_sheet(std::move(sheet));
	}

	IceRegridder *sheet(std::string const &name)
		{ return sheets[sheets_index.at(name)].get(); }

	void filter_cellsA(std::function<bool(long)> const &keepA);

	void filter_cellsA(ibmisc::Domain<int> const &domainA);

	typedef ibmisc::DerefRandomAccessIter<const IceRegridder, typename SheetsT::const_iterator> const_iterator;

	const_iterator begin() const
		{ return const_iterator(sheets.cbegin()); }
	const_iterator end() const
		{ return const_iterator(sheets.cend()); }

	// -----------------------------------------

	/** Adds the weight (native area) of each cell of the Atmosphere grid (A) */
	void wA(SparseVector &w) const;

	/** @return Number of elevation points for a given grid cell */
	unsigned int nhp(int i1) const { return hpdefs.size(); }
	unsigned int nhp() const { return nhp(-1); }
	unsigned long nA() const { return gridA->ndata(); }
	unsigned long nE() const { return nA() * nhp(-1); }

	void ncio(ibmisc::NcIO &ncio, std::string const &vname);

};	// class GCMRegridder
// ===========================================================
// -----------------------------------------------------------
typedef std::function<std::unique_ptr<WeightedSparse>(bool scale)> RegridFunction;

/** Holds the set of "Ur" (original) matrices produced by an IceRegridder. */
class RegridMatrices {
public:
	std::map<std::string, RegridFunction> regrids;

	RegridMatrices(IceRegridder *sheet);

	/** Retrieves a final regrid matrix. */
	std::unique_ptr<WeightedSparse> regrid(std::string const &spec_name, bool scale) const
		{ return (regrids.at(spec_name))(scale); }

};




}	// namespace
