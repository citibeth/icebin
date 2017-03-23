#include <icebin/smoother.hpp>

namespace icebin {

Smoother::~Smoother() {}

// ---------------------------------------------------------
Smoother::Smoother(std::vector<Smoother::Tuple> &&_tuples,
    std::array<double,3> const &_sigma) :
    nsigma(2.),
    nsigma_squared(nsigma*nsigma),
    sigma(_sigma),
    tuples(std::move(_tuples))
{
    // Create an RTree and insert our tuples into it
    for (auto t(tuples.begin()); t != tuples.end(); ++t) {
        rtree.Insert(&t->centroid[0], &t->centroid[0], &*t);
    }
}
// -----------------------------------------------------------
/** Inner loop for Smoother::matrix() */
bool Smoother::matrix_callback(Smoother::Tuple const *t)
{
    // t0 = point from outer loop
    // t = point from innter loop

    // Compute a scaled distance metric, based on the radius in each direction
    double norm_distance_squared = 0;
    for (int i=0; i<3; ++i) {
        double const d = (t->centroid[i] - t0->centroid[i]) / sigma[i];
        norm_distance_squared += d*d;
    }

    if (norm_distance_squared < nsigma_squared) {
        double const gaussian_ij = std::exp(-.5 * norm_distance_squared);
        double w = gaussian_ij * t->area;
        M_raw.push_back(std::make_pair(t->iX_d, w));
        denom_sum += w;
    }
    return true;
}

void Smoother::matrix(TupleListT<2> &ret)
{
    using namespace std::placeholders;  // for _1, _2, _3...

    RTree::Callback callback(std::bind(&Smoother::matrix_callback, this, _1));
    for (auto _t0(tuples.begin()); _t0 != tuples.end(); ++_t0) {
        t0 = &*_t0;
        M_raw.clear();
        denom_sum = 0;

        // Pair t0 with nearby points
        std::array<double,3> min, max;
        for (int i=0; i<3; ++i) {
            min[i] = t0->centroid[i] - nsigma*sigma[i];
            max[i] = t0->centroid[i] + nsigma*sigma[i];
        }
        rtree.Search(min, max, callback);

        // Add to the final matrix
        double factor = 1. / denom_sum;
        for (auto ii=M_raw.begin(); ii != M_raw.end(); ++ii) {
            ret.add({t0->iX_d, ii->first}, factor * ii->second);
        }
    }
}
// -----------------------------------------------------------
void smoothing_matrix(TupleListT<2> &ret_d,
    Grid const *gridX,
    SparseSetT const &dimX,
    DenseArrayT<1> const *elev_s,
    DenseArrayT<1> const &area_d,
    std::array<double,3> const &sigma)
{
    std::vector<Smoother::Tuple> tuples;
    for (auto cell(gridX->cells.begin()); cell != gridX->cells.end(); ++cell) {
        auto iX_s(cell->index);
        if (!dimX.in_sparse(iX_s)) continue;

        auto iX_d(dimX.to_dense(iX_s));

        double area(area_d(iX_d));
        if (area == 0.) (*icebin_error)(-1,
            "Area of cell %ld must be non-zero\n", iX_s);

        double elev((*elev_s)(iX_s));
        if (std::isnan(elev)) (*icebin_error)(-1,
            "Grid cell %ld cannot be masked out\n", iX_s);

        Point const &centroid(gridX->centroid(*cell));
        tuples.push_back(
            Smoother::Tuple(iX_d,
                centroid,
                elev, area));
    }
    Smoother smoother(std::move(tuples), sigma);
    smoother.matrix(ret_d);
}

}    // namespace
