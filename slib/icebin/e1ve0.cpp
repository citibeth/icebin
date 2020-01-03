#include <tuple>
#include <algorithm>
#include <set>
#include <icebin/e1ve0.hpp>

using namespace ibmisc;
using namespace spsparse;

namespace icebin {
namespace e1ve0 {


// ==================================================================
// --------------- Compute E1vE0 based on correction concept

// If you just forgot about the changing basis for vector x, then what
// error would you introduce (on I grid) from that?  The error would
// be:
//         IvE1*x - IvE0*x = [IvE1 - IvE0] x
// We can correct this by adding back the negative.  When converted to
// E1 grid, we get the correction:
//         [E1vI * (IvE0 - IvE1)] x
// Add that to the original value to get corrected value for x (in E space):
//       x + [E1vI * (IvE0 - IvE1)] x
// Remove the x vector to get:
//       E1vE0 = I + [E1vI * (IvE0 - IvE1)]
//
// This will guarantee that E1vE0 is close to I; and E1vE0==I if E1==E0

/** Sort a vector of Tuples by index, and sum entries with same index.
TODO: Move to tuplelist.hpp */
template<class IndexT, class ValT, int RANK>
void consolidate(std::vector<spsparse::Tuple<IndexT,ValT,RANK>> &M)
{
    // Must be at least 2 elements for this to make sense...
    if (M.size() < 2) return;

    // Sort by iX
    std::sort(M.begin(), M.end());

    // Eliminate duplicates
    size_t j=0;
    for (size_t i=1; i<M.size(); ++i) {
        if (M[j].index() == M[i].index()) {
            M[j].value() += M[i].value();
        } else {
            ++j;
            M[j] = M[i];
        }
    }
    M.resize(j+1);
}


spsparse::TupleList<int,double,2> compute_E1vE0c(
std::vector<std::unique_ptr<ibmisc::linear::Weighted_Eigen>> const &XuE1s,  // sparsified
std::vector<std::unique_ptr<ibmisc::linear::Weighted_Eigen>> const &XuE0s,  // sparsified
unsigned long nE,            // Size of (sparse) E vector space, never changes
std::vector<double> const &areaX)
{
    spsparse::TupleList<int,double,2> E1vE0c;
    blitz::Array<double,1> sE1(nE);
    sE1 = 0;

    // 1. Compute UNSCALED correction matrix, per ice sheet
    // 2. Concatenate per-ice-sheet correction matrix into overall correction matrix
    // correct_unscaled = sum_{ice sheet}[ E1uX * (XvE0 - XvE1) ]
    for (size_t i=0; i<XuE1s.size(); ++i) {
        // -------- 1. Compute UNSCALED correction matrix, per ice sheet
        linear::Weighted_Eigen const *XuE1 = &*XuE1s[i];
        linear::Weighted_Eigen const *XuE0 = &*XuE0s[i];

        blitz::Array<double,1> sXuE0(1. / XuE0->wM);
        blitz::Array<double,1> sXuE1(1. / XuE1->wM);

        // blitz::Array<double,1> XuE1s(1. / XuE1->Mw);
        // auto E1vX(map_eigen_diagonal(XuE1s) * XuE1->M->transpose())
        auto E1uX(XuE1->M->transpose());
        auto XvE1(map_eigen_diagonal(sXuE1) * *XuE1->M);
        auto XvE0(map_eigen_diagonal(sXuE0) * *XuE0->M);

        // -------- 2. Concatenate per-ice-sheet correction matrices
        // Assemble the correction matrix, and sparsify its indexing.
        // EigenSparseMatrixT correct_unscaled(E1uX*XvE0 - E1uX*XvE1);
        EigenSparseMatrixT E1vE0c_local(E1uX*(XvE0 - XvE1));
        spcopy(
            accum::ref(E1vE0c),
            E1vE0c_local);

        for (int i=0; i<XuE1->Mw.extent(0); ++i) sE1(i) += XuE1->Mw(i);
    }

    // Convert merged weights to scale factor
    for (int i=0; i<sE1.extent(0); ++i) sE1(i) = 1. / sE1(i);

    // E1vE0c = sE1 * E1vE0c_unscaled
    for (auto ii(E1vE0c.begin()); ii != E1vE0c.end(); ++ii) {
        long const iE1(ii->index(0));
        ii->value() *= sE1(iE1);
    }

    // Tidy up by consolidating
    consolidate(E1vE0c.tuples);

    return E1vE0c;
}

}}
