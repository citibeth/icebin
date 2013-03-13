#pragma once

#include <boost/function.hpp>

namespace glint2 {

// ------------------------------------------------

// ------------------------------------------------

/** Serves two purposes: (1) Translate global grid cell indices to
local indexing scheme, for a particular GCM MPI node. (2) Tells
whether the grid cell is in the node's domain or halo. */
class GridDomain {
public:

	/** Number of indices used locally (by the GCM) to index a grid cell
	(for ModelE, this is generally 2, cells are inexed by i,j) */
	const int num_local_indices;

	GridDomain(int _num_local_indices) : num_local_indies(_num_local_indices) {}
	virtual ~GridDomain() {}

	/** Convert to local indexing scheme for this MPI node. */
	virtual void global_to_local(int gindex_c, int *lindex) const;

	/** Tells whether a grid cell is in the MPI node's domain.
	@param lindex Result of global_to_local() */
	virtual bool in_domain(int *lindex);

	/** Tells whether a grid cell is in the MPI node's halo (or main domain).
	@param lindex Result of global_to_local() */
	virtual bool in_halo(int *lindex);

	/** Default implementation is OK; or re-implement to avoid
	going through extra virtual function call
	@return The in_halo() function */
	virtual boost::function<bool (int)> get_in_halo() const;

	void global_to_local(
		std::vector<int> const &global,
		std::vector<blitz::Array<double,1>> &olocal);

};
// ------------------------------------------------
class GridDomain_Identity {

	GridDomain_Identity() : GridDomain(1) {}

	bool global_to_local(int gindex_c, int *lindex)
		{ lindex[0] = gindex_c; }

	bool in_domain(int index_c) const
		{ return true; }
	bool in_halo(int index_c) const
		{ return true; }
};

// ==========================================================
#if 0
extern void global_to_local(
	GridDomain const &domain,
	std::vector<int> const &global,
	std::vector<blitz::Array<double,1>> &olocal);
#endif

extern std::unique_ptr<VectorSparseMatrix> filter_matrix(
	GridDomain const &domain1,
	GridDomain const &domain2,
	VectorSparseMatrix const &mat);

}
