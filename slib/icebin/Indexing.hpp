#ifndef ICEBIN_INDEXING
#define ICEBIN_INDEXING

namespace icebin {

class Indexing {
public:

	virtual ~Indexing() {};
	virtual long ij_to_index(int i, int j) const = 0;
	virtual void index_to_ij(long index, int &i, int &j) const = 0;
};

} // Namespace
#endif
