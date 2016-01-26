#pragma once

#include <unordered_set>
#include <memory>
#include <blitz/array.h>

#include <ibmisc/netcdf.hpp>

#include <icebin/Grid.hpp>
#include <icebin/sparse.hpp>

namespace icebin {

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
)


class GCMRegridder;

/** Creates regridding matrices for a single ice sheet. */
class IceRegridder {
public:
	BOOST_ENUM_VALUES(Type, int,
		(L0)	(0)
		(L1)	(1)
	)

protected:
	friend class GCMRegridder;

	Type type;
	GCMRegridder const *gcm;	/// Parent pointer
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
	std::unordered_map<long,double> elevI_hash();

#if 0
	// ------------------------------------------------
	/** Number of dimensions of ice vector space */
	virtual long nI() const = 0;

	/** Number of dimensions of interpolation grid vector space. */
	virtual long nG() const
		{ return nI(); }

	// ------------------------------------------------
#endif

	/** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
	NOTE: wAvAp == sApvA */
	void wAvAp(SparseVector &w);

	/** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
	NOTE: wAvAp == sApvA */
	void wEvEp(SparseVector &w);

	void ncio(ibmisc::NcIO &ncio, std::string const &vname);

};	// class IceRegridder

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

	ibmisc::Indexing<int,long> indexingA;
	ibmisc::Domain<long> domainA;				// What's in our MPI halo?
	bool correctA;		/// Should we correct for projection and geometric error?


	/** Convert between (iA, iHP) <--> (iE) */
	ibmisc::Indexing<long,long> indexingHP;
	/** Position of height points in elevation space (same for all GCM
	grid cells) */
	std::vector<double> hpdefs;	// [nhp]


	typedef std::map<std::string, std::unique_ptr<IceRegridder>> SheetsT;
	SheetsT sheets;

public:
	GCMRegridder() {}

	void clear();
	void init(
		std::unique_ptr<Grid> &&_gridA,
		ibmisc::Indexing<int,long> &&_indexingA,
		ibmisc::Domain<long> &&_domainA,		// Tells us which cells in gridA to keep...
		ibmisc::Indexing<long,long> &&_indexingHP,
		bool _correctA);

	// -----------------------------------------

	void add_sheet(std::unique_ptr<IceRegridder> &&sheet)
	{
		sheet->gcm = this;
		sheets.insert(std::make_pair(sheet->name, std::move(sheet)));
	}

	void filter_cellsA(std::function<bool(long)> const &keepA);


	// typedef ibmisc::DerefSecondIter<std::string, IceRegridder, typename SheetsT::iterator> iterator;
	typedef ibmisc::DerefSecondIter<std::string, const IceRegridder, typename SheetsT::const_iterator> const_iterator;

#if 0
	iterator begin()
		{ return iterator(sheets.begin()); }
	iterator end()
		{ return iterator(sheets.end()); }
	const_iterator cbegin() const
		{ return const_iterator(sheets.cbegin()); }
	const_iterator cend() const
		{ return const_iterator(sheets.cend()); }
#endif
	const_iterator begin() const
		{ return const_iterator(sheets.cbegin()); }
	const_iterator end() const
		{ return const_iterator(sheets.cend()); }

	IceRegridder const &at(std::string const &name)
		{ return *sheets.at(name); }

	// -----------------------------------------

	/** Adds the weight (native area) of each cell of the Atmosphere grid (A) */
	void wA(SparseVector &w) { gridA->native_areas(w); }

	/** @return Number of elevation points for a given grid cell */
	int nhp(int i1) const { return hpdefs.size(); }
	int nA() const { return gridA->ndata(); }
	int nE() const { return nA() * nhp(-1); }

	void ncio(ibmisc::NcIO &ncio, std::string const &vname);
};	// class GCMRegridder

}	// namespace
