#include <icebin/smoother.hpp>
#include <spatialindex/capi/sidx_impl.h>
#include <spatialindex/capi/sidx_api.h>
#include <spatialindex/capi/sidx_config.h>

namespace icebin {

Smoother::~Smoother() {}

std::unique_ptr<Index> new_rtree()
{
    using namespace SpatialIndex;

    // create a property set with default values.
    // see utility.cc for all defaults  http://libspatialindex.github.io/doxygen/Utility_8cc_source.html#l00031
    std::unique_ptr<Tools::PropertySet> ps(GetDefaults());
    Tools::Variant var;

    // set index type to R*-Tree
    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = RT_RTree;
    ps->setProperty("IndexType", var);

    // Set index to store in memory (default is disk)
    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = RT_Memory;
    ps->setProperty("IndexStorageType", var);

    // initalise index
    std::unique_ptr<Index> idx(new Index(*ps));

    // check index is ok
    if (!idx->index().isIndexValid()) (*icebin_error)(-1,
        "Error creating SpatialIndex RTree");

    return idx;
}
// ---------------------------------------------------------
struct GaussianWeightsVisitor : public SpatialIndex::IVisitor
{
    // Constants inherited from outer loop
    double const radius_squared;
    double const two_sigma_squared_inv;

    GaussianWeightsVisitor(double _radius_squared, double _two_sigma_squared_inv) :
        radius_squared(_radius_squared),
        two_sigma_squared_inv(_two_sigma_squared_inv)
    {}

    // Things cleared/set each time around outer loop
    Smoother::Tuple *t0;
    std::vector<std::pair<int,double>> M_raw;
    double mass_sum;

    void set_t0(Smoother::Tuple *_t0) {
        t0 = _t0;
        M_raw.clear();
        mass_sum = 0;
    }


    void visitNode(const SpatialIndex:: INode &in) {}
    void visitData(std::vector<const SpatialIndex::IData *> &v) {}
    void visitData(const SpatialIndex::IData &in);
};

void GaussianWeightsVisitor::visitData(const SpatialIndex::IData &in)
{
    Smoother::Tuple *t;
    in.getData((uint32_t)sizeof(t), (byte **)&t);

    double const dx = t->centroid.x - t0->centroid.x;
    double const dy = t->centroid.y - t0->centroid.y;
    double const distance_squared = dx*dx + dy*dy;
    if (distance_squared < radius_squared) {
        double w0 = std::exp(
            -two_sigma_squared_inv * distance_squared);
        M_raw.push_back(std::make_pair(t->iX_d, w0));
        mass_sum += w0 * t->mass;
    }
}


Smoother::Smoother(std::vector<Smoother::Tuple> &&_tuples) :
    tuples(std::move(_tuples))
{
    using namespace SpatialIndex;

    // Create an RTree and insert our tuples into it
    rtree = new_rtree();
    std::array<double,2> coords;
    for (auto t(tuples.begin()); t != tuples.end(); ++t) {
        coords = {t->centroid.x, t->centroid.y};
        Smoother::Tuple *t_ptr = &*t;
        rtree->index().insertData(
            sizeof(t_ptr), &t_ptr,
            SpatialIndex::Point(*(coords[0]),2),
            t->iX_d);
    }
}

void Smoother::matrix(MakeDenseEigenT::AccumT &ret, double sigma)
{
    double radius = 3.*sigma;
    double radius_squared = radius*radius;
    double two_sigma_squared_inv = 1./(2.*sigma*sigma);

    GaussianWeightsVisitor gweights(radius_squared, two_sigma_squared_inv);
    for (auto t0(tuples.begin()); t0 != tuples.end(); ++t0) {
        gweights.set_t0(&*t0);    // Clear accumulators

        // Pair t0 with nearby points
        std::array<double,2> const plow
            {t0->centroid.x - radius, t0->centroid.y - radius};
        std::array<double,2> const phigh
            {t0->centroid.x + radius, t0->centroid.y + radius};
        rtree->intersectsWithQuery(
            SpatialIndex::Region(&plow[0], &phigh[0], 2), gweights);

        // Add to the final matrix
        double factor = 1. / gweights.mass_sum;
        for (auto ii=gweights.M_raw.begin(); ii != gweights.M_raw.end(); ++ii) {
            ret.add({ii->first, t0->iX_d}, factor * ii->second);
        }
    }
}


smoothing_matrix(MakeDenseEigenT::AccumT &ret,
    Grid const &gridX, DenseArrayT<1> const &mass_d, double sigma)
{
    std::vector<Smoother::Tuple> tuples;
    for (int iX_d=0; iX<mass_d.size(); ++iX_d) {
        double mass(mass_d(iX_d));
        if (mass == 0.) continue;

        tuples.push_back(Smoother::Tuple(iX_d,
            gridX.centroid(gridX.at(dimX.to_sparse(iX_d))),
            mass);
    }
    Smoother smoother(std::move(tuples), sigma);
    smoother.matrix(ret);
}

}    // namespace
