#include <icebin/smoother.hpp>

namespace icebin {

Smoother::~Smoother() {}

// ---------------------------------------------------------
Smoother::Smoother(std::vector<Smoother::Tuple> &&_tuples, double sigma) :
    radius(2.*sigma),
    radius_squared(radius*radius),
    two_sigma_squared_inv(1./(2.*sigma*sigma)),
    tuples(std::move(_tuples))
{
    // Create an RTree and insert our tuples into it
    std::array<double,2> coords;
//    std::array<double,2> min,max;
    for (auto t(tuples.begin()); t != tuples.end(); ++t) {
        coords = {t->centroid.x, t->centroid.y};
        rtree.Insert(&coords[0], &coords[0], &*t);
//        double eps=.01;
//        min = {t->centroid.x-eps, t->centroid.y-eps};
//        max = {t->centroid.x+eps, t->centroid.y+eps};
//        rtree.Insert(&min[0], &max[0], &*t);
    }
}
// -----------------------------------------------------------
/** Inner loop for Smoother::matrix() */
bool Smoother::matrix_callback(Smoother::Tuple const *t)
{
    double const dx = t->centroid.x - t0->centroid.x;
    double const dy = t->centroid.y - t0->centroid.y;
    double const distance_squared = dx*dx + dy*dy;
    if (distance_squared < radius_squared) {
        double const gaussian_ij = std::exp(
            -two_sigma_squared_inv * distance_squared);
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
        rtree.Search(
            {t0->centroid.x - radius, t0->centroid.y - radius},
            {t0->centroid.x + radius, t0->centroid.y + radius},
            callback);

        // Add to the final matrix
        double factor = 1. / denom_sum;
        for (auto ii=M_raw.begin(); ii != M_raw.end(); ++ii) {
            ret.add({t0->iX_d, ii->first}, factor * ii->second);
        }
    }
}
// -----------------------------------------------------------
void smoothing_matrix(TupleListT<2> &ret,
    Grid const &gridX, SparseSetT const &dimX,
    DenseArrayT<1> const &area_d, double sigma)
{
    std::vector<Smoother::Tuple> tuples;
    for (int iX_d=0; iX_d<area_d.size(); ++iX_d) {
        double area(area_d(iX_d));
        if (area == 0.) continue;

        tuples.push_back(
            Smoother::Tuple(iX_d,
                gridX.centroid(gridX.cells.at(dimX.to_sparse(iX_d))),
                area));
    }
    Smoother smoother(std::move(tuples), sigma);
    smoother.matrix(ret);
}

}    // namespace
