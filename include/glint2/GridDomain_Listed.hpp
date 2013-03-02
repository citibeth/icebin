#pragma once

#include <glint2/GridDomain.hpp>

namespace glint2 {

// ------------------------------------------------
class GridDomain {

	std::unordered_set<int> domain;
	std::unordered_set<int> halo;

public:

	void halo_add(int index_c)
	{
		halo.add(index_c);
		domain.add(index_c);	// Things in domain are also in halo
	}
	void domain_add(int index_c)
		{ domain.add(index_c); }


	bool in_domain(int index_c) const
		{ return domain.find(index_c) != domain.end(); }

	bool in_halo(int index_c) const
		{ return halo.find(index_c) != halo.end(); }

	// Default implementation is OK; or re-implement to avoid
	// going through extra virtual function call
	boost::function<bool (int)> get_in_halo();
};
// ------------------------------------------------


}
