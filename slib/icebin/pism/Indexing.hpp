#ifndef ICEBIN_PISM_INDEXING_HPP
#define ICEBIN_PISM_INDEXING_HPP

#include <icebin/Indexing.hpp>

namespace icebin {
namespace pism {


class Indexing : public icebin::Indexing {

public:

	const int nx, ny;

	Indexing(int _nx, int _ny) : nx(_nx), ny(_ny) {}


	long ij_to_index(int i, int j) const
		{ return i*ny + j; }	// PISM ordering

	void index_to_ij(long index, int &i, int &j) const
	{
		i = index / ny;
		j = index - i*ny;
	}
};


}}

#endif	// Guard
