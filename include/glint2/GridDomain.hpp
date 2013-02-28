#pragma once

namespace glint2 {

// ------------------------------------------------

// ------------------------------------------------
class GridDomain {
public:
	// # of indices used locally (by the GCM) to index a grid cell
	// (for ModelE, this is generally 2, cells are inexed by i,j)
	const int num_local_indices;

	GridDomain(int _num_local_indices) : num_local_indies(_num_local_indices) {}
	virtual ~GridDomain() {}

	virtual bool global_to_local(int gindex_c, int *lindex) const;
	virtual bool in_domain(int *lindex);
	virtual bool in_halo(int *lindex);

//	virtual bool in_domain(int index_c) const = 0;
//	virtual bool in_halo(int index_c) const = 0;

	// Default implementation is OK; or re-implement to avoid
	// going through extra virtual function call
	virtual boost::function<bool (int)> get_in_halo() const;
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



}
