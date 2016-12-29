
// ------------------------------------------------------
inline std::vector<DenseVector> make_dense_vectors(
    size_t nvars, size_t n)
{
    std::vector<DenseVector> ret;
    for (int i=0; i < nvars; ++i)
        ret.push_back(DenseVector(n));
    return ret;
}

inline std::vector<SparseVector> make_sparse_vectors(
    size_t nvars, size_t n)
{
    std::vector<SparseVector> ret;
    for (int i=0; i < nvars; ++i)
        ret.push_back(SparseVector({n}));
    return ret;
}
// ------------------------------------------------------

