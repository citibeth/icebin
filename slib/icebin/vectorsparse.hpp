#ifndef ICEBIN_VECTORSPARSE_HPP
#define ICEBIN_VECTORSPARSE_HPP

// TODO: Expand and put in ibmisc/linear.  See tuplelist.hpp

namespace icebin {

/** Vector-based sparse matrix with sparse indexing */
template<class IndexT, class ValT, int RANK>
class VectorSparse {
public:
    typedef IndexT index_type;
    typedef ValT val_type;
    static const int rank = RANK;

    // https://stackoverflow.com/questions/4353203/thou-shalt-not-inherit-from-stdvector
    std::vector<IndexT> indices;
    std::vector<ValT> values;

    void clear()
    {
        indices.clear();
        values.clear();
    }

    void add(std::array<index_type,rank> const &index, ValT const &value)
    {
        for (auto i=0; i<rank; ++i) indices.push_back(index[rank]);
        values.push_back(value);
    }

    template<class ArchiveT>
    void serialize(ArchiveT &ar, const unsigned int file_version)
    {
        ar & indices;
        ar & values;
    }

#if 0   // Not a very useful function
    /** Export this sparse matrix to Fortran-supplied variable pointers
        See:
            https://stackoverflow.com/questions/30152073/how-to-pass-c-pointer-to-fortran
    */
    void export_matrix(
        int *&indices_p,
        double *&values_p,
        int &nele)
    {
        nele = values.size();
        indices_p = &indices[0];
        values_p = &values[0];
    }
#endif

};

} // namespace icebin
#endif
