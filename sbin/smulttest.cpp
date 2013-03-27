#include <iostream>
#include <blitz/array.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <giss/SparseMatrix.hpp>

blitz::Array<double,2> random_matrix(int m, int n)
{
	static boost::random::mt19937 gen;

	// Make one matrix
    boost::random::uniform_int_distribution<> idist(0, 2);
    boost::random::uniform_real_distribution<> rdist(-50000.0,10000.0);

	blitz::Array<double,2> aa(m, n);
	for (int i=0; i<aa.extent(0); ++i) {
	for (int j=0; j<aa.extent(1); ++j) {
		int id = idist(gen);
		if (id != 0) aa(i,j) = 0.0;
		else aa(i,j) = rdist(gen);
//printf("[%d, %d] ==> %g\n", i, j, aa(i,j));
	}}

	return aa;
}

blitz::Array<double, 2> multiply(blitz::Array<double, 2> const &aa, blitz::Array<double, 2> const &bb)
{
	// aa_ij, bb_jk
	int ni = aa.extent(0);
	int nk = bb.extent(1);
	int nj = aa.extent(1);
	if (nj != bb.extent(0)) throw std::exception();

	blitz::Array<double,2> cc(ni,nk);
	for (int i=0; i<ni; ++i) {
	for (int k=0; k<nk; ++k) {
		double cval = 0;
		for (int j=0; j<nj; ++j) cval += aa(i,j) * bb(j,k);
		cc(i,k) = cval;
	}}

	return cc;
}

giss::VectorSparseMatrix dense_to_sparse(blitz::Array<double,2> aa, int ncopy)
{
	giss::VectorSparseMatrix ret(giss::SparseDescr(
		aa.extent(0), aa.extent(1)));

	double ncopy_inv = 1.0 / (double)ncopy;
	for (int i=0; i<aa.extent(0); ++i) {
	for (int k=0; k<aa.extent(1); ++k) {
		if (aa(i,k) == 0) continue;
		double val = aa(i,k) * ncopy_inv;
		for (int x=0; x<ncopy; ++x) ret.add(i, k, val);
	}}

	return ret;
}

void print(giss::VectorSparseMatrix const &cc_s)
{
	printf("-------------- BEGIN Sparse\n");
	for (auto ii=cc_s.begin(); ii != cc_s.end(); ++ii) {
		printf("[%d %d] --> %g\n", ii.row(), ii.col(), ii.val());
	}
	printf("-------------- END Sparse\n");
}

int main(int argc, char **argv)
{
	int size[3] = {3,4,5};

	auto aa(random_matrix(size[0], size[1]));
	auto bb(random_matrix(size[1], size[2]));

	auto cc(multiply(aa,bb));


	std::cout << aa << std::endl;
	std::cout << bb << std::endl;
	std::cout << cc << std::endl;


	auto aa_s(dense_to_sparse(aa,1));
	auto bb_s(dense_to_sparse(bb,2));

//	print(aa_s);
//	print(bb_s);

	auto cc_sp(giss::multiply(aa_s, bb_s));
	print(*cc_sp);
}
