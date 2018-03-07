#ifndef ICEBIN_RLWEIGHTEDSPARSE_HPP
#define ICEBIN_RLWEIGHTEDSPARSE_HPP

#include <blitz/array.h>
#include <ibmisc/zsparsearray.hpp>
#include <ibmisc/netcdf.hpp>

namespace icebin {

BOOST_ENUM_VALUES(SparseFillType, int,
    (zero) (0)
    (nan) (1)
)

struct RLWeightedSparse {
    // All sparse indexing here!

    ibmisc::ZSparseArray<int,double,1> wM;            // ==0 in M's nullspace
    ibmisc::ZSparseArray<int,double,2> M;    // BvA
    ibmisc::ZSparseArray<int,double,1> Mw;            // ==0 in M's nullspace

    /** True if this regridding matrix is conservative.  Matrices could be
    non-conservative, for example, in the face of smoothing on I.  Or when
    regridding between the IceBin and ModelE ice sheets. */
    bool conservative = false;


    RLWeightedSparse(std::array<long,2> shape)
        : wM({shape[0]}), M(shape), Mw({shape[1]}) {}

    /** Applies a regrid matrix.
    ALL INDEXING IS SPARSE!!!

    Nominally computes B{in} = BvA{ij} * A{jn}
    (where BvA is this)
    If the matrix ix non-conservative, it also accounts for that.

        |i| = Size of B vector space
        |j| = Size of A vector space
        |n| = Number of vectors being transformed

    @param A The values to regrid, as a series of row vectors
    */
    void apply(
        // this = BvA
        blitz::Array<double,2> const &A,      //  IN: A{nj} one row per variable
        blitz::Array<double,2> &B,            // OUT: B{ni} one row per variable
        SparseFillType fill_type=SparseFillType::zero,
        bool force_conservation=true) const;

    void ncio(ibmisc::NcIO &ncio, std::string const &vname);

};


};    // namespace icebin
#endif    // guard
