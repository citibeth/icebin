#ifndef ICEBIN_E1VE0_HPP
#define ICEBIN_E1VE0_HPP

#include <icebin/eigen_types.hpp>


/** Constructs the E1vE0 matrix, used to regrid elevation class-based
state variables (eg land surface state) from one coupling timestep to
the next. */

namespace icebin {
namespace e1ve0 {

// One basis function: Tuple of just (iI,value) part of EvI;
// represents just the basis function
typedef std::vector<spsparse::Tuple<long,double,1>> BasisFn;

/** Set of basis functions, arranged by gridcell and elevation class. */
typedef std::map<
    long /*iE*/,
    std::pair<long /*iA*/, BasisFn>
> BasisFnMap;

// Convenience accessors...
inline long const &iE(BasisFnMap::const_iterator const &ii)
    { return ii->first; }
inline long const &iA(BasisFnMap::const_iterator const &ii)
    { return ii->second.first; }
inline BasisFn const &bfn(BasisFnMap::const_iterator const &ii)
    { return ii->second.second; }
inline BasisFn &bfn(BasisFnMap::iterator const &ii)
    { return ii->second.second; }



/** Computes innerproduct area basis functions, A gridcell by A gridcell.
@param E0vIs Original E0vIs matrices, converted to basis function form
@param areaX Area of (concatenated) exchange grid
@return List of ovlerap of each pair: (iE0, iE0, value) */
extern spsparse::TupleList<long,double,2> compute_E1vE0_scaled(
BasisFunctionMap const &bfn1s,
BasisFunctionMap const &bfn0s,
unsigned long nE,            // Size of (sparse) E vector space, never changes
std::vector<double> const &areaX)



/** Extracts basis functions from a set of EvX matrices.
@param XvEs XvE matrix for each ice sheet (X = exchange grid)
@param basesI Offset to add to iX (local) to convert to iX (global).
*/
extern BasisFnMap extract_basis_fns(
std::vector<SparseSetT const *> const &dimEs,
std::vector<EigenSparseMatrixT const *> const &XvEs,
std::vector<long> const &basesX)


}}    // namespace
#endif    // guard
