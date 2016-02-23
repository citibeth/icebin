




template<class SparseMatrixT, class EigenSparseMatrixT>
void eigen_to_matrix(
	SparseMatrixT &ret,
	EigenSparseMatrixT const &M,
	std::array<SparseSet<
		typename SparseMatrixT::index_type,
		typename EigenTripletT::Index> *, 2> const &dims)
{
	for (int k=0; k<M.outerSize(); ++k) {
	for (EigenSparseMatrixT::InnerIterator ii(M,k); ii; ++ii) {
		ret.add({dims[0]->to_sparse(ii.row(), dims[1]->to_sparse(ii.col())}, ii.value());
	}}
}





template <class SparseMatrixT>
class SparseTriplets
{
public:
	typedef Eigen::Triplet<typename SparseMatrixT::val_type> EigenTripletT;
	typedef Egien::SparseMatrix<typename SparseMatrixT::val_type> EigenSparseMatrixT;

	SparseMatrixT M;
	std::array<SparseSet<
		typename SparseMatrixT::index_type,
		typename EigenTripletT::Index> *, 2> const dims;

	void add(
		std::array<typename SparseMatrixT::index_type, 2> index,
		typename SparseMatrixT::val_type val)
	{
		M.push_back(index[0], index[1], val);
		dims[0]->add(index[0]);
		dims[1]->add(index[1]);
	}

	EigenSparseMatrixT to_eigen(char transpose) const
	{
		std::vector<EigenTripletT> triplets;

		for (auto ii=M.begin(); ii != M.end(); ++ii) {
			auto dense0 = dims[0]->to_dense(transpose == 'T' ? 1 : 0);
			auto dense1 = dims[1]->to_dense(transpose == 'T' ? 0 : 1);
			triplets.push_back(EigenTripletT(dense0, dense1, ii.val()));
		}
		EigenSparseMatrixT ret(
			dims[transpose == 'T' ? 1 : 0].size(),
			dims[transpose == 'T' ? 0 : 1].size());
		ret.setFromTriplets(triplets.begin(), triplets.end());
		return ret;
	}
};
