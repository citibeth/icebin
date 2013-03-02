#include <glint2/GridDomain.hpp>

namespace glint2 {

// -------------------------------------------------------
boost::function<bool (int)> GridDomain::get_in_halo()
{ return boost::bind(this, &GridDomain::in_halo, _1); }

void global_to_local(
	GridDomain const &domain,
	std::vector<int> const &global,
	std::vector<blitz::Array<double,1>> &olocal)
{
	if (olocal.size() != domain.num_local_indices) {
		fprintf(stderr, "MatrixDomainer::get_rows() had bad dimension 1 = %d (expected %ld)\n", olocal.extent(1), domain.num_local_indices);
		throw std::exception();
	}

	for (auto ii = olocal.begin(); ii != olocal.end(); ++ii) {
		// Make sure it has the right dimensions
		if (olocal[i].extent(0) != global.size()) {
			fprintf(stderr, "MatrixDomainer::get_rows() had bad dimension 0 = %d (expected %ld)\n", olocal.extent(0), global.size());
			throw std::exception();
		}
	}


	// Copy out data, translating to local coordinates
	for (int j=0; j < global.size(); ++j) {
		int lindex[domain.num_local_indices];
		domain.global_to_local(global[j], lindex);
		for (int i=0; i<domain.num_local_indices; ++i)
			olocal[i](j) = lindex[i];
	}
}
// -------------------------------------------------------
std::unique_ptr<VectorSparseMatrix> filter_matrix(
	GridDomain const &domain1,
	GridDomain const &domain2,
	VectorSparseMatrix const &mat)
{
	VectorSparseMatrix ret;

	for (auto ii = mat.begin(); ii != mat.end(); ++ii) {
		// Output of linear transformation: Only include
		// if it's part of our domain.
		if (!domain1.in_domain(ii.row())) continue;

		// Input of linear transformation: must be in halo
		if (!domain2.in_halo(ii.col())) {
			fprintf("Error filtering matrix: grid cell %d in input (column) is not available in the halo.\n", ii.col());
			throw std::exception;
		}

		ret.add(ii.row(), ii.col(), ii.val());
	}

	return ret;
}
// -------------------------------------------------------


}
