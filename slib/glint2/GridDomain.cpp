#include <glint2/GridDomain.hpp>

namespace glint2 {

// -------------------------------------------------------
boost::function<bool (int)> GridDomain::get_in_halo2() const
	{ return boost::bind(&GridDomain::in_halo2, this, _1); }

#if 0
void GridDomain::global_to_local(
	blitz::Array<double,1> const &global,
	std::vector<blitz::Array<double,1>> &olocal)
{
	if (olocal.size() != this->num_local_indices) {
		fprintf(stderr, "MatrixDomainer::get_rows() had bad dimension 1 = %d (expected %ld)\n", olocal.extent(1), this->num_local_indices);
		throw std::exception();
	}

	for (auto ii = olocal.begin(); ii != olocal.end(); ++ii) {
		// Make sure it has the right dimensions
		if (olocal[i].extent(0) != global.extent(0)) {
			fprintf(stderr, "MatrixDomainer::get_rows() had bad dimension 0 = %d (expected %ld)\n", olocal.extent(0), global.extent(0));
			throw std::exception();
		}
	}


	// Copy out data, translating to local coordinates
	for (int j=0; j < global.extent(0); ++j) {
		int lindex[this->num_local_indices];
		this->global_to_local(global(j), lindex);
		for (int i=0; i<this->num_local_indices; ++i)
			olocal[i](j) = lindex[i];
	}
}
#endif
// -------------------------------------------------------
std::unique_ptr<giss::VectorSparseMatrix> filter_matrix(
	GridDomain const &domain1,
	GridDomain const &domain2,
	giss::VectorSparseMatrix const &mat)
{
	std::unique_ptr<giss::VectorSparseMatrix> ret(
		new giss::VectorSparseMatrix(mat));

	int lindex1[domain1.num_local_indices];
	int lindex2[domain2.num_local_indices];
	for (auto ii = mat.begin(); ii != mat.end(); ++ii) {
		// Output of linear transformation: Only include
		// if it's part of our domain.
		domain1.global_to_local(ii.row(), lindex1);
		if (!domain1.in_domain(lindex1)) continue;

		// Input of linear transformation: must be in halo
		domain2.global_to_local(ii.col(), lindex2);
		if (!domain2.in_halo(lindex2)) {
			fprintf(stderr, "Error filtering matrix: grid cell %d in input (column) is not available in the halo.\n", ii.col());
			throw std::exception();
		}

		ret->add(ii.row(), ii.col(), ii.val());
	}

	return ret;
}
// -------------------------------------------------------


}
