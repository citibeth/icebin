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

        /** Position of the "center" of the grid cell.  Assumption is
        grid cells are small, they are weighted according to the value
        of the Gaussian function at this centroid. */
        Point const centroid;

        /** Mass to assign to the grid cell; typically the area /
        integral of this cell, as it overlaps some other grid. */
        double mass;

        Tuple(int _iX_d,
            Point const _centroid, double _mass)
        : iX_d(_iX_d), centroid(_centroid), mass(_mass)
        {}
    };

    typedef ibmisc::RTree<Tuple const *, double, 2> RTree;

    double const radius;
    double const radius_squared;
    double const two_sigma_squared_inv;

protected:
    // ------ Variables used in loop; see Smoother::matrix()
    Tuple *t0;
    std::vector<std::pair<int,double>> M_raw;
    double mass_sum;
    // --------------------------

    RTree rtree;

    // Outer loop: set once
    std::vector<Tuple> tuples;

public:
    /** Set up the RTree needed for the smoothing matrix */
    Smoother(std::vector<Tuple> &&_tuples, double sigma);

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
    @param ret
        Receptacle for the smoothing matrix.
    @param gridX
        Grid being used (sparse indexing)
    @param dimX
        Translation between dense and sparse indexing
    @param mass_d
        Mass (area for L0 grid) of each grid cell.  This
        will typically be a weight vector for a regridding matrix; for
        example, for IvA or IvE.  There is mass only for grid cells
        that overlap the A/E grid.
    @param sigma
        Scale size of the Gaussian smoothing function.  If smoothing
        the IvA matrix, than sigma should be about the size of an A
        grid cell.
*/
extern void smoothing_matrix(TupleListT<2> &ret,
    Grid const &gridX, SparseSetT const &dimX,
    DenseArrayT<1> const &mass_d, double sigma);

}

#endif
