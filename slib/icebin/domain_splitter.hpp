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

// ----------------------------------------------------------------------
/** Scales and merges (unscaled) matrices stored in WTs */
TupleListLT<2> merge_scale(std::vector<WeightedTuple const *> const &WTs)
{
    auto shape(WTs[0].shape());

    // Create weight array by merging weights in the tuples
    blitz::Array<double,1> wM(shape[0]);
        wM = 0;
    for (TupeListLT<2> const *Mi : M0) {
    for (auto &tp : Mi->wM.tuples) {
        wM(tp.index(0)) += tp.value();
    }}

    // Convert weights to scale factors
    for (int i=0; i<wM.extent(0); ++i) {
        wM(i) = 1. / wM(i);    // Will produce Inf where unused
    }

    // Rename wM to sM to reflect change of meaning
    blitz::Array<double,1> sM;
    sM.reference(wM);

    // Copy the main matrix, scaling!
    TupeListLT<2> M;    // Output
    M.set_shape(shape);
    for (TupeListLT<2> const *Mi : M0) {
    for (auto &tp : Mi->tuples) {
        M.add(tp.index(), tp.value() * sM(tp.index(0)));
    }}

    return M;
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
    for (int iAE=0; iAE != (int)IndexAE::COUNT; ++iAE) {
        auto &gcm_ivalsX(out.gcm_ivalss_s[iAE]);
        DomainDecomposerT const &domainsX(*domainsAE[iAE]);

        for (size_t i=0; i<out.gcm_ivalss_s[iAE].index.size(); ++i) {
            // Convert from index to tuple
            long const iA(gcm_ivalsX.index[i]);

            // Figure out MPI domain of the resulting tuple
            int idomain = domainsX.get_domain(iA);

            // Add our values to the appropriate MPI domain
//printf("domain1 outs %ld %ld %d %d\n", iA, outs.size(), idomain, nvar[iAE]);
//printf("domain1 %d %ld", iAE, iA, 
            outs[idomain].gcm_ivalss_s[iAE].add(iA, &gcm_ivalsX.vals[i*nvar[iAE]]);
        }
    }

    // Split the standard sparse vectors and matrices
    // Copying a TupleListLT<2>
    spcopy(
        accum::domain(strideb, &domainsE,
        &outs.front().E1vE0_s),
        out.E1vE0_s)

    return outs;
}

}    // namespace icebin
#endif
