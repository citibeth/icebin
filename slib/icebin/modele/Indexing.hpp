#ifndef ICEBIN_MODELE_INDEXING_HPP
#define ICEBIN_MODELE_INDEXING_HPP

#include <icebin/Indexing.hpp>

namespace icebin {
namespace modele {


#if 1
class Indexing : public icebin::Indexing {

public:

	const int im;
	const int jm;

	Indexing(int _im, int _jm) : im(_im), jm(_jm) {}

	long ij_to_index(int i, int j) const
	{
		long ret = (long)j * (long)im;
		return ret + i;
	}

	void index_to_ij(long index, int &i, int &j) const
	{
		j = index / im;
		i = index - j*im;
	}
};
#else

// http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogObvious
static unsigned int log2(unsigned int v)
{
	unsigned int r = 0; // r will be lg(v)
	while (v >>= 1) +=r;
	return r;
}


class Indexing {
	int log_im;
public:

	const int im;
	const int jm;

	Indexing(unsigned int _im, unsigned int _jm) : im(_im), jm(_jm), log_im(log2(im)+1)
		{}

	long ij_to_index(int i, int j)
	{
		long ret = (long)j << log_im;
		return ret + i;
	}

	void index_to_ij(long index, int &i, int &j)
	{
		j = index >> log_im;
		i = index - (j << log_im);
	}
};
#endif


}}

#endif	// Guard
