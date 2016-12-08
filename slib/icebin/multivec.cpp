#include <icebin/multivec.hpp>
#include <icebin/error.hpp>

namespace icebin {

template<class IndexT, class ValT>
void VectorMultivec::add(IndexT const &ix, ValT const *val)
{
    index.push_back(ix);
    for (int i=0; i<nvar; ++i) vals.push_back(val[i]);
}

VectorMultivec concatenate(std::vector<VectorMultivec> const &vecs)
{
    VectorMultivec ret(vecs[0].size());

    if (vecs.size() == 0) (*icebin_error)(-1,
        "Must concatenate at least one vector");

    ret.nvar = vecs[0].nvar;

    for (auto ii=vecs.begin(); ii != vecs.end(); ++ii) {
        ret.index.insert(ret.index.end(), ii->index.begin(), ii->index.end());
        ret.vals.insert(ret.vals.end(), ii->vals.begin(), ii->vals.end());

        if (ret.nvar != ii->nvar) (*icebin_error)(-1,
            "Inconsistant nvar: %d vs %d", ret.nvar, ii->nvar);
    }
    return ret;
}

}
