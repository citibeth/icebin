#ifndef ICEBIN_SMOOTHER_HPP
#define ICEBIN_SMOOTHER_HPP

#include <ibmisc/RTree.hpp>
#include <icebin/IceRegridder.hpp>

namespace icebin {

/** Producer of smoothing matrices that is (somewhat) independent of
    particular grid implementations.  Smoother::Tuple encapsulates all
    we need to know about each grid cell to make the smoothing
    matrix. */
class Smoother {
public:
    /** Required information about each grid cell */
    struct Tuple {
        /** A unique index for this grid cell. */
        int iX_d;    // Dense index for this cell

        /** Position (xyz) of the "center" of the grid cell.  Assumption is
        grid cells are small, they are weighted according to the value
        of the Gaussian function at this centroid. */
        std::array<double,3> centroid;

        /** Area to assign to the grid cell; typically the area /
        integral of this cell, as it overlaps some other grid. */
        double area;

        Tuple(int _iX_d,
            std::array<double,2> _centroid, double z, double _area)
        : iX_d(_iX_d),
        centroid{_centroid[0], _centroid[1], z},
        area(_area)
        {}
    };

    typedef ibmisc::RTree<Tuple const *, double, 3> RTree;

    // Include points out to nsigma*sigma away
    std::array<double,3> sigma;
    double const nsigma;
    double const nsigma_squared;

protected:
    // ------ Variables used in loop; see Smoother::matrix()
    Tuple *t0;
    std::vector<std::pair<int,double>> M_raw;
    double denom_sum;
    // --------------------------

    RTree rtree;

    // Outer loop: set once
    std::vector<Tuple> tuples;

public:
    /** Set up the RTree needed for the smoothing matrix */
    Smoother(std::vector<Tuple> &&_tuples, std::array<double,3> const &sigma);

    ~Smoother();    // Not inline because of forward-declared type of rtree

protected:
    /** Inner loop for Smoother::matrix() */
    bool matrix_callback(Tuple const *t);

public:
    /** Generate the smoothing matrix. */
    void matrix(TupleListT<2> &ret);
};

/** Produces a smoothing matrix that "smears" one grid cell into
    neighboring grid cells, weighted by a radial Gaussian.
    The vector space (gridX) in question has two forms of indexing:
        * Sparse (_s): Native indexing, no guarantee that indices are closely spaced.
          Some grid cells could be unused.
        * Dense (_d): A dense renumbering of grid cells, resulting in ONLY the ones needed
          for this (and other current) operations.
    @param ret_d
        Receptacle for the smoothing matrix.
        Uses dense indexing.
    @param gridX
        Grid being used (sparse indexing)
    @param dimX
        Translation between dense (local to our computing) and sparse
        (grid-native) indexing.
    @param elev_s
        Elevation of each grid cell.
        NaN for grid cells that are masked out.
        Indexed by sparse (grid-native) indexing
    @param area_d
        (Generalized) area of each grid cell.  This
        will typically be a weight vector for a regridding matrix; for
        example, for IvA or IvE.  There is area only for grid cells
        that overlap the A/E grid.
        Indexed by dense indexing
    @param sigma
        Scale size of the Gaussian smoothing function.  If smoothing
        the IvA matrix, than sigma[0] and sigma[1] should be about the
        size of an A grid cell and sigma[2] infinity.  If smoothing IvE,
        then sigma[2] should be about the elevation difference between
        different elevation classes.
*/
extern void smoothing_matrix(TupleListT<2> &ret,
    AbbrGrid const &agridX,
    SparseSetT const &dimX,
    DenseArrayT<1> const &elev_s,
    DenseArrayT<1> const &area_d,
    std::array<double,3> const &sigma);

}

#endif
