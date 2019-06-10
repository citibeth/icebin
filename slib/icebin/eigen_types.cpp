#include <spsparse/SparseSet.hpp>
#include <icebin/eigen_types.hpp>
#include <ibmisc/linear/compressed.hpp>

using namespace spsparse;

namespace icebin {

EigenSparseMatrixT to_eigen_M(  // (generates in dense indexing)
ibmisc::ZArray<int,double,2> const &BvA,
std::array<SparseSetT *,2> dims)
{
    // ======================= Create a merged EOpvAOp of base ice and ice sheets
    // (and then call through to _compute_AAmvEAm)

    // Accumulator for merged M (unscaled)
    MakeDenseEigenT BvA_m(
        {SparsifyTransform::ADD_DENSE},   // convert sparse to dense indexing
        dims, '.');
    auto BvA_accum(BvA_m.accum());

    // Copy elements to accumulator matrix
    for (auto ii(BvA.generator()); ++ii; ) {
        BvA_accum.add({ii->index(0), ii->index(1)}, ii->value());
    }

    // Set overall size
    std::array<long,2> BvA_shape(BvA.shape());
    dims[0]->set_sparse_extent(BvA_shape[0]);
    dims[1]->set_sparse_extent(BvA_shape[1]);

    // Compute M and wAOp
    return BvA_m.to_eigen();
}

}    // namespace
