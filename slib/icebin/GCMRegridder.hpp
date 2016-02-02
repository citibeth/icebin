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
	int nhp(int i1) const { return hpdefs.size(); }
	int nA() const { return gridA->ndata(); }
	int nE() const { return nA() * nhp(-1); }

	void ncio(ibmisc::NcIO &ncio, std::string const &vname);
};	// class GCMRegridder
// ===========================================================
template<class TypeT>
class Transpose {
public:
	TypeT const * const M;
	bool const transpose;

	Transpose(TypeT const *_M, bool _transpose) : M(_M), transpose(_transpose) {}
#if 0
	Transpose(Transpose<TypeT> &&rhs) : M(rhs.M), transpose(rhs.transpose) {}
	void operator=(Transpose<TypeT> &&rhs) {
		M = rhs.M;
		transpose = rhs.transpose;
	}
#endif
};
// -----------------------------------------------------------
/** Helper class... */
template<class TypeT>
class VectorAllocator {
	std::vector<std::unique_ptr<TypeT>> mem;
public:
	TypeT *alloc() {
		std::unique_ptr<TypeT> M(new TypeT());
		TypeT *ret = M.get();
		mem.push_back(std::move(M));
		return ret;
	}
};
// -----------------------------------------------------------
/** Holds the set of "Ur" (original) matrices produced by an IceRegridder. */
class RegridMatrices {
protected:
	IceRegridder *sheet;

	VectorAllocator<SparseMatrix> matrix_mem;
	VectorAllocator<SparseVector> vector_mem;

	std::map<std::string, Transpose<SparseMatrix>> regrids;
	std::map<std::string, SparseVector const *> diags;

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
	void regrid(SparseMatrix &ret, std::string const &spec_name) const;

	/** Retrieves a weight or scale vector.  Options are:
		w/sAvG		Partial cell weighting for A
		w/sEvG		Partial cell weighting for E
		w/sA		Whole cell weighting for A
		w/sG		Whole cell weighting for E
	*/
	void weight(SparseVector &ret, std::string const &spec_name) const;

};




}	// namespace
