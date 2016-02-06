#pragma once

#include <unordered_set>
#include <memory>
#include <blitz/array.h>

#include <ibmisc/netcdf.hpp>

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
	std::unordered_map<long,double> elevI_hash();

	// ------------------------------------------------
	/** Number of dimensions of ice vector space */
	virtual long nI() const = 0;

	/** Number of dimensions of interpolation grid vector space. */
	virtual long nG() const = 0;

	// ------------------------------------------------

	/** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
	NOTE: wAvAp == sApvA */
	void wAvAp(SparseVector &w);

	/** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
	NOTE: wAvAp == sApvA */
	void wEvEp(SparseVector &w);

	virtual void GvEp_noweight(
		SparseMatrix &ret,
		std::unordered_map<long,double> const &elevIh) const = 0;

	virtual void GvI_noweight(SparseMatrix &ret, std::unordered_map<long,double> const &elevIh) const = 0;
	virtual void GvAp_noweight(SparseMatrix &ret) = 0;

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

	ibmisc::Domain<int> domainA;				// What's in our MPI halo?
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
		ibmisc::Domain<int> &&_domainA,		// Tells us which cells in gridA to keep...
		std::vector<double> &&_hpdefs,
		ibmisc::Indexing<long,long> &&_indexingHP,
		bool _correctA);

	// -----------------------------------------

	void add_sheet(std::unique_ptr<IceRegridder> &&sheet)
	{
		sheet->gcm = this;
		sheets.insert(std::make_pair(sheet->name, std::move(sheet)));
	}

	IceRegridder *sheet(std::string const &name)
		{ return sheets.at(name).get(); }

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
	void wA(SparseVector &w) const;

	/** @return Number of elevation points for a given grid cell */
	unsigned int nhp(int i1) const { return hpdefs.size(); }
	unsigned long nA() const { return gridA->ndata(); }
	unsigned long nE() const { return nA() * nhp(-1); }

	void ncio(ibmisc::NcIO &ncio, std::string const &vname);
};	// class GCMRegridder
// ===========================================================
template<class TypeT>
class Transpose {
public:
	TypeT M;
	bool const transpose;

	Transpose(TypeT const &&_M, bool _transpose) : M(std::move(_M)), transpose(_transpose) {}
};
// -----------------------------------------------------------
template<class TypeT>
class LazyPtr {
	TypeT *_ptr;						// Borrowed reference
	std::unique_ptr<TypeT> _uptr;	// Owned reference
	std::function<std::unique_ptr<TypeT> ()> _compute_owned;
	std::function<TypeT *()> _compute_borrowed;
public:
	LazyPtr() {}

	LazyPtr(std::unique_ptr<TypeT> &&ptr) : _uptr(std::move(ptr)) {}

	LazyPtr(std::function<std::unique_ptr<TypeT> ()> &&compute_owned)
		: _compute_owned(std::move(compute_owned)) {}

	LazyPtr(std::function<TypeT *()> &&compute_borrowed)
		: _compute_borrowed(std::move(compute_borrowed)) {}

	LazyPtr(TypeT *ptr) : _ptr(ptr) {}

	TypeT &operator*() const {
		if (!_ptr) {
			if (_compute_borrowed) {
				const_cast<TypeT *>(_ptr) = _compute_borrowed();
			} else {
				const_cast<std::unique_ptr<TypeT>>(_uptr) = _compute_owned();
				const_cast<TypeT *>(_ptr) = &*_uptr;
			}
		}
		return *ptr;
	}
	TypeT *operator->() const
		{ return &operator*(); }
};
// -----------------------------------------------------------
/** Holds the set of "Ur" (original) matrices produced by an IceRegridder. */
class RegridMatrices {
protected:
	IceRegridder *sheet;

	std::map<std::string, Transpose<LazyPtr<SparseMatrix>>> regrids;
	std::map<std::string, LazyPtr<SparseVector>> diags;

	SparseVector *invert(SparseVector *v);

	/** Adds a regrid matrix and its variants.
	@param G The destination vector space of the matrix.  This must always be "G".
	@param Z The source vectors space of the matrix.
	@param m_GvZ Memory where to store the underlying matrix. */
	void add_regrid(
		std::string const &G,
		std::string const &Z,
		SparseMatrix *GvZ);


	void add_weight(
		std::string const &B,
		std::string const &A,
		std::function<void(SparseVector &)> const &scale_fn);


	/** Definitions of the Ur regridding and scaling matrices that
	must be multiplied together to obtain the final regridding
	matrices. */
	static std::map<std::string, std::array<std::string, 5>> regrid_specs;

	// ----------------------------------------------------------
	void wAvG(SparseVector &ret) const;
	void wEvG(SparseVector &ret) const;
	void wA(SparseVector &ret) const
		{ sheet->gcm->wA(ret); }

	// using namespace std::placeholders;
	static std::map<std::string, std::function<void(RegridMatrices const *, SparseVector &)>> weight_specs;
	// ----------------------------------------------------------
		

public:

	RegridMatrices(IceRegridder *sheet);

	/** Retrieves a final regrid matrix. */
	void regrid(SparseMatrix &ret, std::string const &spec_name, Weighting weighting) const;

	/** Retrieves a weight or scale vector.  Options are:
		w/sAvG		Partial cell weighting for A
		w/sEvG		Partial cell weighting for E
		w/sA		Whole cell weighting for A
		w/sG		Whole cell weighting for E
	*/
	void weight(SparseVector &ret, std::string const &spec_name) const;

};




}	// namespace
