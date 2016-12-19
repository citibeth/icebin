#ifndef ICEBIN_DOMAIN_SPLITTER_HPP
#define ICEBIN_DOMAIN_SPLITTER_HPP

#include <icebin/GCMCoupler.hpp>

namespace icebin {

// ----------------------------------------------------------------
// SAMPLE to go in DomainDecomposerT
// class DomainDecomposer {
// public:
// 
//     /** Number of domains */
//     size_t size();
// 
//     /** Returns the MPI rank of grid cell */
//     int get_rank(long ix);
// }
// ----------------------------------------------------------------
// split_by_domain() implementation...
//
/** Accumulator converts 1-D indices to n-D indices */
template<class AccumulatorT, class DomainDecomposerT>
class DomainAccum {
    std::vector<AccumulatorT *> subs;    // One per domain
    DomainDecomposerT *domains;

public:
    SPSPARSE_LOCAL_TYPES(AccumulatorT);    // Must be rank=1

    DomainAccum(
        std::vector<AccumulatorT *> &&_subs,
            DomainDecomposerT *_domains)
        : subs(std::move(_subs)), domains(_domains)
    {
        if (rank != 1) (*icebin_error)(-1,
            "DomainAccum must be instantiated with AccumulatorT of rank=1");
    }

    void set_shape(std::array<long, AccumulatorT::rank> const &shape)
    {
        for (size_t i=0; i<subs.size(); ++i) subs[i]->set_shape(shape);
    }

    void add(std::array<long, AccumulatorT::rank> iA, val_type const &val)
    {
        int iDomain = domains->get_domain(iA[0]);
        subs[iDomain]->add(iA, val);
    }
};
// ----------------------------------------------------------------
/** @param sub0 Fist sub-accumulator in array of accumulators
@param strideb Stride (in bytes) between first and next accumulator in non-dense. */
template<class AccumulatorT, class DomainDecomposerT>
DomainAccum<AccumulatorT,DomainDecomposerT> domain_accum(
    AccumulatorT *sub0, size_t strideb,
    DomainDecomposerT *domains)
{
    std::vector<AccumulatorT *> subs;
    for (size_t i=0; i < domains->size(); ++i) {
        subs.push_back(reinterpret_cast<AccumulatorT *>(
            reinterpret_cast<char *>(sub0) + strideb));
    }
    return DomainAccum<AccumulatorT,DomainDecomposerT>(
        std::move(subs), domains);
}
// ----------------------------------------------------------------------
template<class DomainDecomposerT>
std::vector<GCMInput> split_by_domain(
    GCMInput const &out,
    DomainDecomposerT const &domainsA,
    DomainDecomposerT const &domainsE);

template<class DomainDecomposerT>
std::vector<GCMInput> split_by_domain(
    GCMInput const &out,
    DomainDecomposerT const &domainsA,
    DomainDecomposerT const &domainsE)
{
    // Put domain decomposers in a nice array
    std::array<DomainDecomposerT const *, GridAE::count> domainsAE
        {&domainsA, &domainsE};
    int ndomains = domainsA.size();

    // Construct the output
    std::array<int, GridAE::count> nvar(out.nvar());

    std::vector<GCMInput> outs;
    outs.reserve(ndomains);
    for (size_t i=0; i<ndomains; ++i)
        outs.push_back(GCMInput(nvar));

    size_t strideb = sizeof(GCMInput);

    // Split each element of parallel sparse vectors
    for (int iAE=0; iAE != GridAE::count; ++iAE) {
        auto &gcm_ivalsX(out.gcm_ivalsAE[iAE]);
        DomainDecomposerT const &domainsX(*domainsAE[iAE]);

        for (size_t i=0; i<out.gcm_ivalsAE[iAE].index.size(); ++i) {
            // Convert from index to tuple
            long const iA(gcm_ivalsX.index[i]);

            // Figure out MPI domain of the resulting tuple
            int idomain = domainsX.get_domain(iA);

            // Add our values to the appropriate MPI domain
            outs[idomain].gcm_ivalsAE[iAE].add(iA, &gcm_ivalsX.vals[i*nvar[iAE]]);
        }
    }

    // Split the standard sparse vectors and matrices
    auto E1vE0(domain_accum(&outs.front().E1vE0, strideb, &domainsE));
        copy(E1vE0, out.E1vE0);

    auto AvE1(domain_accum(&outs.front().AvE1, strideb, &domainsA));
        copy(AvE1, out.AvE1);

    auto wAvE1(domain_accum(&outs.front().wAvE1, strideb, &domainsA));
        copy(wAvE1, out.wAvE1);

    auto elevE1(domain_accum(&outs.front().elevE1, strideb, &domainsA));
        copy(elevE1, out.elevE1);

    return outs;
}

}    // namespace icebin
#endif
