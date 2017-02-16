#ifndef ICEBIN_MULTIVEC_HPP
#define ICEBIN_MULTIVEC_HPP

#include <cstddef>
#include <vector>

namespace boost {
namespace serialization {
    class access;
}}

namespace icebin {

/** Stores multiple sparse vectors that use the same grid cells. */
class VectorMultivec {
    friend class boost::serialization::access;
public:
    typedef long index_type;
    typedef std::vector<double> val_type;
    static const int rank = 1;

    // Stores a bunch of parallel sparse vectors
    // Index of each element in the parallel vectors
    std::vector<long> index;
    // Values of for each element in the vectors.  
    std::vector<double> vals;
    // Number of _vals element per _ix element
    int nvar;

    // Needed by boost::mpi::gather() [GCMCoupler_ModelE.cpp]
    VectorMultivec() : nvar(-1) {}

    VectorMultivec(int _nvar) : nvar(_nvar) {}

    /** Number of elements in each parallel array. */
    inline size_t size() const { return index.size(); }

    template<typename ArchiveT>
    void serialize(ArchiveT& ar, const unsigned version) {
        ar & index;
        ar & vals;
        ar & nvar;
    }

    /** Adds a new element to all the sparse vectors */
    void add(long ix, double const *val);

    void add(long ix, std::vector<double> &val)
        { add(ix, &val[0]); }

    double val(int varix, long ix) const
        { return vals[ix*nvar + varix]; }

};

VectorMultivec concatenate(std::vector<VectorMultivec> const &vecs);


}    // namespace icebin

#endif
