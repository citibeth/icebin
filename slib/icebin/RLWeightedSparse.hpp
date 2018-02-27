#ifndef ICEBIN_RLWEIGHTEDSPARSE_HPP
#define ICEBIN_RLWEIGHTEDSPARSE_HPP

#incluce <blitz/array.h>

namespace ibmisc {
    template<class CountT, class IndexT, class ValueT, int RANK>
    class RLSparseArray;
}

namespace icebin {

BOOST_ENUM(SparseFillType
    (ZERO) (0)
    (NAN) (1)
)

struct RLWeightedSparse {
    // All sparse indexing here!

    RLSparseArray<int,int,double,1> const &wM,            // ==0 in M's nullspace
    RLSparseArray<int,int,double,2> const &M,    // BvA
    RLSparseArray<int,int,double,1> const &Mw,            // ==0 in M's nullspace

    std::array<int,2> shape {0,0};

    /** True if this regridding matrix is conservative.  Matrices could be
    non-conservative, for example, in the face of smoothing on I.  Or when
    regridding between the IceBin and ModelE ice sheets. */
    bool conservative = false;


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
    RLWeightedSparse::apply(
        // this = BvA
        blitz::Array<double,2> const &A,      //  IN: A{nj} one row per variable
        blitz::Array<double,2> &B,            // OUT: B{ni} one row per variable
        SparseFillType fill_type=SparseFillType::ZERO,
        bool force_conservation=true);

    void ncio(NcIO &ncio, std::string const &vname)
    {
        auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
        get_or_put_att(info_v, ncio.rw, "conservative", conservative);

        wM.ncio(vname + ".wM", "int", "int", "double");
        M.ncio(vname + ".M", "int", "int", "double");
        mW.ncio(vname + ".Mw", "int", "int", "double");
    }

};


};    // namespace icebin
#endif    // guard
