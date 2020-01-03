#ifndef ICEBIN_E1VE0_HPP
#define ICEBIN_E1VE0_HPP

#include <icebin/eigen_types.hpp>
#include <ibmisc/linear/eigen.hpp>

/** Constructs the E1vE0 matrix, used to regrid elevation class-based
state variables (eg land surface state) from one coupling timestep to
the next. */

namespace icebin {
namespace e1ve0 {

/** Computes the E1vE0c correction matrix to regrid snow/firn model from
one set of elevation classes to the next.

NOTE: Actual matrix applied is (E1vE0 - I), based on correction factor idea
      E1vE0 = <Identity Matrix> + E1vE0c.
      Phrasing it in this way makes it easier a matrix of the type (I + C)
@param XuE1s Latest set of per-ice-sheet XuE matrices (unscaled) (X = exchange grid)
@param XuE0s Previous coupling-timestep set of XuE matrices
@param nE Number of theoretical elevation classes (in sparse E indexing)
@param areaX Area of each exchange gridcell. */
extern TupleListT<2> compute_E1vE0c(
std::vector<std::unique_ptr<ibmisc::linear::Weighted_Eigen>> const &XuE1s,
std::vector<std::unique_ptr<ibmisc::linear::Weighted_Eigen>> const &XuE0s,
unsigned long nE,            // Size of (sparse) E vector space, never changes
std::vector<double> const &areaX);


}}    // namespace
#endif    // guard
