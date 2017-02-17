#ifndef ICEBIN_DOMAIN_SPLITTER_HPP
#define ICEBIN_DOMAIN_SPLITTER_HPP

#include <icebin/GCMCoupler.hpp>
#include <spsparse/accum.hpp>

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

namespace accum {


/* This accumulator splits between domains based on the first dimension (0) of AccumulatorT. */
template<class AccumulatorT, class DomainDecomposerT>
class Domain {
    std::vector<AccumulatorT *> subs;    // One per domain
    DomainDecomposerT *domain_splitter;

public:
    SPSPARSE_LOCAL_TYPES(AccumulatorT);    // Must be rank=1

    Domain(
        std::vector<AccumulatorT *> &&_subs,
            DomainDecomposerT *_domain_splitter)
        : subs(std::move(_subs)), domain_splitter(_domain_splitter)
    {
#if 0
        if (rank != 1) (*icebin_error)(-1,
            "accum::Domain must be instantiated with AccumulatorT of rank=1 (it has rank %d instead)", (int)rank);
#endif

    }

    void set_shape(std::array<long, AccumulatorT::rank> const &shape)
    {
        for (size_t i=0; i<subs.size(); ++i) subs[i]->set_shape(shape);
    }

    void add(std::array<long, AccumulatorT::rank> iA, val_type const &val)
    {
        int iDomain = domain_splitter->get_domain(iA[0]);
        subs[iDomain]->add(iA, val);
    }
};
// ----------------------------------------------------------------

/** @param sub0 Fist sub-accumulator in array of accumulators
@param strideb Stride (in bytes) between first and next accumulator in non-dense. */
template<class AccumulatorT, class DomainDecomposerT>
Domain<AccumulatorT,DomainDecomposerT> domain(
    size_t strideb,
    DomainDecomposerT *domains,
    AccumulatorT *sub0)
{
    std::vector<AccumulatorT *> subs;
    for (size_t i=0; i < domains->size(); ++i) {
        subs.push_back(reinterpret_cast<AccumulatorT *>(
            reinterpret_cast<char *>(sub0) + strideb));
    }
printf("domain() template created %ld sub-domains\n", subs.size());
    return Domain<AccumulatorT,DomainDecomposerT>(
        std::move(subs), domains);
}

}    // namespace accum
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
    using namespace spsparse;

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
        auto &gcm_ivalsX(out.gcm_ivalsAE_s[iAE]);
        DomainDecomposerT const &domainsX(*domainsAE[iAE]);

        for (size_t i=0; i<out.gcm_ivalsAE_s[iAE].index.size(); ++i) {
            // Convert from index to tuple
            long const iA(gcm_ivalsX.index[i]);

            // Figure out MPI domain of the resulting tuple
            int idomain = domainsX.get_domain(iA);

            // Add our values to the appropriate MPI domain
//printf("domain1 outs %ld %ld %d %d\n", iA, outs.size(), idomain, nvar[iAE]);
//printf("domain1 %d %ld", iAE, iA, 
            outs[idomain].gcm_ivalsAE_s[iAE].add(iA, &gcm_ivalsX.vals[i*nvar[iAE]]);
        }
    }

    // Split the standard sparse vectors and matrices
    // Copying a TupleListLT<2>
#   define SPCOPY(MATRIX, DOMAINSX) \
    spcopy( \
        accum::domain(strideb, &DOMAINSX, \
        &outs.front().MATRIX), \
        out.MATRIX)

    SPCOPY(E1vE0_s, domainsE);
    SPCOPY(AvE1_s, domainsA);
    SPCOPY(wAvE1_s, domainsA);
    SPCOPY(elevE1_s, domainsE);

#   undef SPCOPY

    return outs;
}

}    // namespace icebin
#endif
