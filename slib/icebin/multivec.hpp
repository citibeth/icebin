#ifndef ICEBIN_MULTIVEC_HPP
#define ICEBIN_MULTIVEC_HPP

#include <cstddef>
#include <vector>
#include <ibmisc/blitz.hpp>

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
    // Weight (eg grid cell size) corresponding to each index
    std::vector<double> weights;
    // Values of for each element in the vectors.  
    std::vector<double> vals;
    // Number of _vals element per _ix element
    int nvar;

    // Needed by boost::mpi::gather() [GCMCoupler_ModelE.cpp]
    VectorMultivec() :  nvar(-1) {}

    VectorMultivec(int _nvar) : nvar(_nvar) {}

    /** Number of elements in each parallel array. */
    inline size_t size() const { return index.size(); }

    template<typename ArchiveT>
    void serialize(ArchiveT& ar, const unsigned version) {
        ar & index;
        ar & weights;
        ar & vals;
        ar & nvar;
    }

    /** Adds a new element to all the sparse vectors */
    void add(long ix, double const *val, double weight);

    void add(long ix, std::vector<double> &val, double weight)
        { add(ix, &val[0], weight); }

    double val(int varix, long ix) const
        { return vals[ix*nvar + varix]; }


    void to_dense_scale(blitz::Array<double,1> &scaleE) const;

    /** @param scaleE Multiply by this.
    @param denseE Pre-allocated array to put it in */
    void to_dense(
        int ivar,
        blitz::Array<double,1> const &scaleE,
        double fill,
        blitz::Array<double,1> &denseE) const;

};

VectorMultivec concatenate(std::vector<VectorMultivec> const &vecs);


}    // namespace icebin

#endif
