#ifndef ICEBIN_SEARISE_INDEXING_HPP
#define ICEBIN_SEARISE_INDEXING_HPP

#include <icebin/Indexing.hpp>

namespace icebin {
namespace searise {


class Indexing : public icebin::Indexing {

public:

	const int nx, ny;

	Indexing(int _nx, int _ny) : nx(_nx), ny(_ny) {}


	long ij_to_index(int i, int j) const
		{ return j*nx + i; }	// SEARISE ordering

	void index_to_ij(long index, int &i, int &j) const
	{
		j = index / nx;
		i = index - j*nx;
	}
};


}}

#endif	// Guard
