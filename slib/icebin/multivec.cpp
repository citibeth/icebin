#include <icebin/multivec.hpp>
#include <icebin/error.hpp>

namespace icebin {

static double const nan = std::numeric_limits<double>::quiet_NaN();

void VectorMultivec::add(long ix, double const *val, double weight)
{
    index.push_back(ix);
    weights.push_back(weight);
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
        ret.weights.insert(ret.weights.end(), ii->weights.begin(), ii->weights.end());
        ret.vals.insert(ret.vals.end(), ii->vals.begin(), ii->vals.end());

        if (ret.nvar != ii->nvar) (*icebin_error)(-1,
            "Inconsistant nvar: %d vs %d", ret.nvar, ii->nvar);
    }
    return ret;
}

void VectorMultivec::to_dense_scale(blitz::Array<double,1> &scaleE) const
{
   int nE(scaleE.extent(0));

    // Fill our dense scale
    scaleE = 0;
    for (unsigned int i=0; i<this->index.size(); ++i) {
        auto iE(this->index[i]);
        if (iE >= nE) (*icebin_error)(-1,
            "Index out of range: %ld vs. %ld", (long)iE, (long)nE);
        scaleE(iE) += this->weights[i];
    }

    for (int iE=0; iE<nE; ++iE)
        scaleE(iE) = 1. / scaleE(iE);
}


/** @param scaleE Multiply by this.
@param denseE Pre-allocated array to put it in */
void VectorMultivec::to_dense(
    int ivar,
    blitz::Array<double,1> const &scaleE,
    double fill,
    blitz::Array<double,1> &denseE) const
{
    int nE(denseE.extent(0));

    // Fill our dense var
    denseE = nan;   // nan = This has not been touched yet
    for (unsigned int i=0; i<this->index.size(); ++i) {
        auto iE(this->index[i]);
        if (iE >= nE) (*icebin_error)(-1,
            "Index out of range: %ld vs. %ld", (long)iE, (long)nE);
        if (std::isnan(denseE(iE))) {
            denseE(iE) = this->val(ivar, i) * scaleE(iE);
        } else {
            denseE(iE) += this->val(ivar, i) * scaleE(iE);
        }
    }

    // Set all untouched items to the fill value
    for (int iE=0; iE<denseE.extent(0); ++iE) {
        if (std::isnan(denseE(iE))) denseE(iE) = fill;
    }

}

}
