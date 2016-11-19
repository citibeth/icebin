
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
#if 0
        // Get the CSR sparse matrix to convert GCM outputs to ice model inputs
        ibmisc::CSRAndUnits trans = vt.apply_scalars({
            std::make_pair("by_dt", 1.0 / ((itime - api->itime_last) * api->dtsrc)),
            std::make_pair("unit", 1.0)});
#endif


SparseSetT dimA, dimE;
dimA.add_sorted(AvE.dim_begin(0), AvE.dim_end(0));
dimE.add_sorted(AvE.dim_begin(1), AvE.dim_end(1));
SparseMatrix AdvE;    // Ad = Dense A grid...


/** Algorithm remaps one dimension of a VectorCooArray, to densify the
    indices.  Usage:

    VectorCooArray<long, double, 2> IvE(...);
    SparseSet<long,long> dimE;
    dimE.add_sorted(IvE.dim_begin(1), IvE.dim_end(1));
    SparseMatrix IvEd;    // Ed = Densified indices on E grid...
    densify_one_dim(IvEd, IvE, dimE, 1);

        or alternately in-place:

    OverwriteAccum<decltype(IvE)::iterator> IvEd(IvE.begin());
    densify_one_dim(IvEd, IvE, dimE, 1);
*/
template<class VectorCooArrayT, class AccumulatorT, class SparseSetT>
void densify_one_dim(AccumulatorT &ret, VectorCooArrayT const &A, SparseSetT const &map, int map_dim)
{
    std::array<int,VectorCooArrayT::rank> idx;
    for (auto ii=A.begin(); ii != A.end(); ++ii) {
        VectorCooArrayT::indices_type index(ii.index);
        index[map_dim] = map.to_dense(index[map_dim]);
        ret.add(index, ii.val());
    }
}


// std::vector<SparseVector> transform_and_regrid(
// SparseMatrix const &M,     // Regrids
// CSRMatrix const &T,        // Transorms logical variables
// std::vector<blitz::Array<double,1>> const &inn,
// std::vector<blitz::Array<double,1>> const &out)
// {
//     const size_t nvar_inn = inn.size();
//     const size_t nvar_out = out.size();
//     blitz::Array<double,1> zval(nvar_out);    // Transformed value of inn
// 
//     // Do the multiplication
//     for (auto iM = M.begin(); iM != M.end(); ++iM) {
// 
//         // Transform units on the input while multiplying by M
//         zval = 0;
//         for (int xi=0; xi<nvar_out; ++xi) {
//             std::vector<std::pair<int, double>> const &row(trans.mat[xi]);
//             for (auto xjj=row.begin(); xjj != row.end(); ++xjj) {
//                 int xj = xjj->first;
//                 double xval = xjj->second;
//                 out[xi](iM.index(0)) += iM.val() * xval * inn[xj](iM.index(1));
//             }
//         }
//     }
//         
// }
// 
// 
