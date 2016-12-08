#include <icebin/GCMCoupler.hpp>



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
    SPSPARSE_LOCAL_TYPES(AccumulatorT);    // Must be rank=1

    std::vector<AccumulatorT *> subs;    // One per domain
    DomainDecomposerT *domains;

public:

    DomainAccum(
        std::vector<AccumulatorT *> &&_subs,
        ibmisc::Indexing *_indexing,
            DomainDecomposerT *_domains)
        : subs(std::move(_subs)), indexing(_indexing), domains(_domains) {}

    void set_shape(std::array<long, rank> const &shape)
    {
        for (size_t i=0; i<subs.size(); ++i) subs[i]->set_shape(shape);
    }

    void add(long iA, value_type const &val)
    {
        int iDomain = domain->get_domain(iA);
        subs[iDomain]->add(iA, val);
    }
};
// ----------------------------------------------------------------
/** @param sub0 Fist sub-accumulator in array of accumulators
@param strideb Stride (in bytes) between first and next accumulator in non-dense. */
template<class AccumulatorT, class DomainDecomposerT>
void domain_accum(AccumulatorT *sub0, size_t strideb,
        DomainDecomposerT *domains)
{
    std::vector<AccumulatorT *> subs;
    for (size_t i=0; i < domains->size(); ++i) {
        subs.push_back(reinterpret_cast<AccumulatorT *>(
            reinterpret_cast<char *>(sub0) + strideb);
    }
    return DomainTupleAccum<AccumulatorT,DomainDecomposerT>(
        std::move(subs), indexing, domains);
}
// ----------------------------------------------------------------------
template<class DomainDecomposerT>
std::vector<GCMCouplerOutput> split_by_domain(
    GCMCouplerOutput const &out,
    DomainDecomposerT const &domainsA,
    DomainDecomposerT const &domainsE);

template<class DomainDecomposerT>
std::vector<GCMCouplerOutput> split_by_domain(
    GCMCouplerOutput const &out,
    DomainDecomposerT const &domainsA,
    DomainDecomposerT const &domainsE)
{
    // Put domain decomposers in a nice array
    std::array<DomainDecomposer<IndexT> *, GridAE::count> domainsAE
        {&domainsA, domainsE};
    int ndomains = domainsA.size();
    int nvar = out.nvar();

    // Construct the output
    std::vector<GCMCouplerOutputT> outs;
    outs.reserve(ndomains);
    for (size_t i=0; i<ndomains; ++i)
        outs.push_back(GCMCouplerOutputT(nvar));

    size_t strideb = sizeof(GCMCouplerOutput);

    // Split each element of parallel sparse vectors
    for (GridAE iAE=0; iAE != GridAE::count; ++iAE) {
        auto &gcm_ivalsX(gcm_ivals[iAE]);
        DomainDecomposer &domainsX(domainsAE[iAE]);

        for (size_t i=0; i<out.gcm_ivals[iAE].index.size(); ++i) {
            // Convert from index to tuple
            long const iA(gcm_ivalsX.index[i]);

            // Figure out MPI domain of the resulting tuple
            int idomain = domainX.domain(iA);

            // Add our values to the appropriate MPI domain
            outs[idomain].add(iA, &gcm_ivalsX.vals[i*nvar]);
        }
    }

    // Split the standard sparse vectors and matrices
    auto E1vE0(domain_accum(&outs.front().E1vE0, strideb, &domainsE));
        copy(E1vE0, out.E1vE0);

    auto AvE1(domain_tuple_accum(&outs.front().AvE1, strideb, &domainsA));
        copy(AvE1, out.AvE1);

    auto wAvE1(domain_tuple_accum(&outs.front().wAvE1, strideb, &domainsA));
        copy(wAvE1, out.wAvE1);

    auto elevE1(domain_tuple_accum(&outs.front().elevE1, strideb, &domainsA));
        copy(elevE1, out.elevE1);

    return outs;
}

// ------------------------------------------------------------
