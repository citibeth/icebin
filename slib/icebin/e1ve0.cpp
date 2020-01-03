#include <tuple>
#include <algorithm>
#include <set>
#include <icebin/e1ve0.hpp>

namespace icebin {
namespace e1ve0 {

/** Computes innerproduct of two basis functions, as extracted from the matrix.
This is done based on the innerproduct of their common gridcells in I.
@param bfn0 List of (iI,value) of basis function.  Sorted and consolidated.
@param bfn1 List of (iI,value) of basis function.  Sorted and consolidated.
@param areaX Area of (concatenated) exchange grid */
static double innerproduct0(
BasisFn const &bfn0,    // std::vector<spsparse::Tuple<long,double,1>>
BasisFn const &bfn1,    // std::vector<spsparse::Tuple<long,double,1>>
std::vector<double> const &areaX)
{
    double innerproduct = 0;

    auto ii0(bfn0.begin());
    auto ii1(bfn1.begin());

    while (true) {
        while (ii0->index(0) < ii1->index(0)) {
            ++ii0; if (ii0 == bfn0.end()) goto end_loop;
        }
        while (ii1->index(0) < ii0->index(0)) {
            ++ii1; if (ii1 == bfn1.end()) goto end_loop;
        }

        // Elements where we match!
        // Assumes: only one element of each index
        while (ii0->index(0) == ii1->index(0)) {
            long const iX = ii0->index(0);

            innerproduct += areaX[iX] * ii0->value() * ii1->value();

            ++ii0; if (ii0 == bfn0.end()) goto end_loop;
            ++ii1; if (ii1 == bfn1.end()) goto end_loop;
        }
    }

end_loop:
    return innerproduct;
}


/** Computes innerproduct area of basis functions within a single A gridcell.
It does this by trying the Cartesian products of all basis functions.
@param areaX Area of (concatenated) exchange grid
@param innerproducts OUT: Store inner products between basis functions here. */
static void innerproductEs(
BasisFunctionMap::iterator const &ii0a,
BasisFunctionMap::iterator const &ii0b,
blitz::Array<double,1> &wE0,
BasisFunctionMap::iterator const &ii1a,
BasisFunctionMap::iterator const &ii1b,
std::vector<double> const &areaX,
TupleListT<2> &innerproducts)
{
    for (auto ii0(ii0a); ii0!=ii0b; ++ii0) {
        for (auto ii1(ii1a); ii1!=ii1b; ++ii1) {
            double const val = innerproduct0(bfn(ii0), bfn(ii1), areaX);
            if (val > 0) {
                innerproducts.add(
                    std::array<long,2>{iE(ii0)), iE(ii1)}, val);
                wE0(iE(ii0)) += val;
            }
        }
    }
}

/** Computes innerproduct area basis functions, A gridcell by A gridcell.
@param E0vIs Original E0vIs matrices, converted to basis function form
@param areaX Area of (concatenated) exchange grid
@return List of ovlerap of each pair: (iE0, iE0, value) */
spsparse::TupleList<long,double,2> compute_E1vE0_scaled(
BasisFunctionMap const &bfn1s,
BasisFunctionMap const &bfn0s,
unsigned long nE,            // Size of (sparse) E vector space, never changes
std::vector<double> const &areaX)
{
    TupleListT E1vE0(std::array<long,2>{nE,nE});

    // Initialize weights array
    blitz::Array<double,0> wE1(nE);
    wE1 = 0;

    // Innerproduct sets of basis functions, by each A gridcell
    BasisFnMap::iterator ii1b(bfn1s.begin());
    for (auto ii1a(ii1b); ii1b != bfn1s.end(); ii1a=ii1b) {
        // Find end of the current chunk of (iA,*)
        // NOTE: ii->second.first is iA
        for (auto ii1b=ii1a+1;
            ii1b != bfn1s.end() && iA(ii1b) == iA(ii1a); ++ii1b);

        BasisFnMap::iterator ii0b(bfn0s.begin());

        // Do Cartesian product of all basis functions in this array
        for (auto ii0a(ii0b); ii0b != bfn0s.end(); ii0a=ii0b) {
            // Find end of the current chunk of (iA,*)
            for (auto ii0b=ii0a+1;
                ii0b != bfn0s.end() && iA(ii0b) == iA(ii0a); ++ii0b);

            // Do Cartesian inner products
            innerproductEs(ii1a,ii1b,wE1,  ii0a,ii0b, areaX, E1vE0);
        }
    }

    // Convert weights to scale, and apply to result
    for (size_t i=0; i<nE; ++i) wE1(i) = 1. / wE1(i);
    for (auto ii=E1vE0.begin(); ii != E1vE0.end(); ++ii)
        ii->val() *= wE1(ii->index(0));

    return E1vE0;
}

/** Extracts basis functions from a set of EvX matrices.
@param XvEs XvE matrix for each ice sheet (X = exchange grid)
@param basesI Offset to add to iX (local) to convert to iX (global).
*/
BasisFnMap extract_basis_fns(
std::vector<SparseSetT const *> const &dimEs,
std::vector<EigenSparseMatrixT const *> const &XvEs,
std::vector<long> const &basesX)
{
    BasisFnMap ret;

    // Combine all matrices into one set of basis functions...
    for (size_t iM=0; iM < XvEs.size(); ++iM) {
        // Dig through the matrix, separating it into basis functions.
        long last_iE = -1;    // TODO: This caching will only be effective if XvEs are sorted in that way
        BasisFnMap::iterator iiE;
        for (ii=begin(*XvEs[iM]); ii != ent(*XvEs[iM]); ++ii) {
            long const iX = ii->index(0) + basesX[iM];
            long const iEd = ii->index(1);
            auto iE = dimEs[iM].to_sparse(iEd);

            // Look up basis function (iiE) corresponding to iE
            if (iE != last_iE) {
                auto const tup(indexingHC.index_to_tuple(iE));
                    long const iA = tup[0];
                    long const ihc = tup[1];

                iiE = ret.find(iE);
                if (iiE == ret.end()) {
                    auto xx(ret.insert(std::make_pair(iA, BasisFn())));
                    iiE = xx.first;    // insert() returns iterator to inserted item
                }
                last_iE = iE;
            }

            // Add element to the basis function
            bfn(iiE).push_back(spsparse::Tuple<long,double,1>(std::array<long,1>{iX}, ii->value()));
        }
    }
}

}}
