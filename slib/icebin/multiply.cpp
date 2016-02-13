#include <icebin.sparse.hpp>
#include <spsparse/SparseSet.hpp>
#include <Eigen/SparseCore>

using namespace spsparse;

namespace icebin {

typedef SparseSet<long,int> IBSparseSet;
typedef Eigen::SparseVector<double> ESparseVector;
typedef Eigen::SparseMatrix<double> ESparseMatrix;
typedef Eigen::Triplet<double> ETriplet;

ESparseMatrix diag_to_eigen(SparseVector const &V, IBSparseSet const &dim)
{
	// Construct our set of triplets
	std::vector<ETriplet> triplets;
	for (auto ii=V.begin(); ii != V.end(); ++ii) {
		auto sparse_index = ii.index(0);
		auto dense_index = dim.to_dense(sparse_index);
		triplets.push_back(ETriplet(dense_index, dense_index, ii.val()));
	}

	// Construct and return a diagonal Eigen sparse matrix from the triplets
	ESparseMatrix ret(dim.size(), dim.size());
	ret.setFromTriplets(triplets.begin(), triplets.end());
	return ret;
}

ESparseMatrix matrix_to_eigen(SparseVector const &V, IBSparseSet const &dim0, IBSparseSet const &dim1)
{
	// Construct our set of triplets
	std::vector<ETriplet> triplets;
	for (auto ii=V.begin(); ii != V.end(); ++ii) {
		auto dense0 = dim0.to_dense(ii.index(0));
		auto dense1 = dim1.to_dense(ii.index(1));
		triplets.push_back(ETriplet(dense0, dense1, ii.val()));
	}

	// Construct and return an Eigen sparse matrix from the triplets
	ESparseMatrix ret(dim0.size(), dim1.size());
	ret.setFromTriplets(triplets.begin(), triplets.end());
	return ret;
}

void eigen_to_matrix(SparseMatrix &ret, ESparseMatrix const &M)
{
	for (int k=0; k<M.outerSize(); ++k) {
	for (SparseMatrix<double>::InnerIterator ii(mat,k); ii; ++ii) {
		ret.add({ii.row(), ii.col()}, ii.value());
	}}
}


/** Computes: diag(scales[0]) * diag(scales[1]) * M */
void multiply(SparseMatrix &ret, std::vector<SparseVector *> const &scales, SparseMatrix const &M)
{

	// Translate to dense indices
	std::array<SparseSet<long,int>, 2> dims;

	// Don't worry abut sorting, just add everything into our dense dimension 0
	for (size_t k=0; k<scales.size(); ++k)
		dims[0].add(scales[k].dim_begin(0), scales[k]->dim_end(0));
	dims[0].add(M.dim_begin(0), M.dim_end(0));
	dims[1].add(M.dim_begin(1), M.dim_end(1));

	// Densify our dimensions, while converting to Eigen::SparseMatrix
	std::vector<ESparseMatrix> emats;
	emats.push_back(diag_to_eigen(*scales[0], dims[0]));
	emats.push_back(diag_to_eigen(*scales[1], dims[0]));
	emats.push_back(matrix_to_eigen(M, dims[0], dims[1]));

	// Multiply the matrices together
	ESparseMatrix eret = emats[0];
	for (int k=1; k<emats.size(); ++k) eret = eret * emats[i];

	// Undensify our result
	eigen_to_matrix(ret, eret);
}

}
