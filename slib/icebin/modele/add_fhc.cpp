#include <spsparse/eigen.hpp>
#include <spsparse/blitz.hpp>


/** Arrays from TOPO (and GIC), BEFORE we started modifying them. */
class TopoFile {
    icebin::Grid *gridA;

    // From TOPO
    blitz::Array<double,1> focean, flake, fgrnd, fgice, zatmo, hlake;
    // From GIC
    blitz::Array<double,1> zatmo_m;
};



/** Sort and remove duplicates in a std::vector */
template<class TypeT, class SumT>
static void consolidate(std::vector<TypeT> &vals, SumT const &addto)
{
    std::sort(vals.begin(); vals.end());
    size_t i=1;
    size_t j=1;
    TypeT last = vals[0];
    for (; i < vals.size(); ++i) {
        if (vals[i] != last) vals[j++] = vals[i];
        else addto(vals[j], vals[i]);
    }
    vals.resize(i);
}






WeightedSparse const &AvE;
WeightedSparse const &EvI_nc;    // No projection correction
WeightedSparse const &AvI;
WeightedSparse const &IvE;

DenseVector<1> const &elevI


// C++ Translation of add_fhc.py
// This is the body of an update_topo() subroutine.

/**
@param elevIs Current elevI for each ice sheet (including NaN for
    masked-out cells).
@param areaA Native area of the grid cells (ModelE-style sparse indexing)
*/
void add_fhc(
    GCMCoupler const *gcm_coupler,
    IceRegridder const *ice_regridder,

    // Inputs
    std::vector<blitz::Array<double,1>> elevIs,
    areaA  // native_area

OUTPUTS:
//    spsparse::TupleList<long,double,1> elevE(nE());
//    spsparse::TupleList<long,double,1> elevA(nA());

    // Eg = GCM's Elevation Space
    // E = IceBin's Elevation Space
    blitz::Array<double,1> elevEg(gcm_coupler->nE_gcm());   // All EC's, not just icebin-related
    blitz::Array<double,1> elevA(gcm_regridder->nA());
{

    std::vector<std::array<long,1>> indicesEg;    // Indices we've seen
    std::vector<std::array<long,1>> indicesA;
    blitz::Array<double,1> wA(nA());

    HCSegmentData const &legacy_segment(gcm_coupler->gcm_params.segment("legacy"));
    HCSegmentData const &sealand_segment(gcm_coupler->gcm_params.segment("sealand"));
    HCSegmentData const &ec_segment(gcm_coupler->gcm_params.segment("ec"));

    auto const nhc_ice(ec_segment.size);
    auto const icebin_base_hc(ec_segment.base);
    long E_to_Eg = ec_segment.base * nA();    // Offset into elevE global array

    wA = 0;
    for (size_t i=0; i < gcm_regridder->ice_regridders.size(); ++i) {
        SparseSetT dimA, dimE;
        TmpAlloc tmp;

        DenseVector<1> const &elevI(*elevIs[i]);
        IceRegridder const *ice_regridder(&*gcm_regridder->ice_regridders[i]);

        ice_regridder->set_elevI(elevI);

        // A SparseSet that is identity for the entire range of I
        SparseSetT dimI;
        auto dimI(id_sparse_set<SparseSetT>(ice_regridder->nI()));

        RegridMatrices::Params rm_params(
            true,        // scale
            true,        // correctA
            {0,0,0},     // sigma (smoothing)
            true);       // conserve
        RegridMatrices rm(ice_regridder);

        auto EvI(rm.matrix("EvI", {&dimE, &dimI}, rm_params));
        auto AvI(rm.matrix("AvI", {&dimA, &dimI}, rm_params));
        auto AvE(rm.matrix("AvE", {&dimA, &dimE}, rm_params));


        // Copy sparse indices into global vectors
        for (auto ii(dimE.sparse_begin()); ii != dimE.sparse_end(); ++ii)
            indicesE.push_back(*ii);
        for (auto ii(dimA.sparse_begin()); ii != dimA.sparse_end(); ++ii)
            indicesA.push_back(*ii);


        // Copy into global wA array
        // (Icebin-covered area of each A grid cell)
        spcopy(
            accum::to_sparse(make_array(&dimA),
            accum::blitz_existing(wA)),
            AvE.wM);    // blitz::Array<double,1>


        // Add to global elevE
        {
            auto elevEi(EvI.apply(elevI, 0));
            spcopy(
                accum::to_sparse(make_array(&dimE),
                accum::add_index({E_to_Eg},
                accum::blitz_existing(elevEg))),
                elevEi, false);
        }


        // Add to global elevA
        {
            auto elevAi(AvI.apply(elevI, 0));
            spcopy(
                accum::to_sparse(make_array(&dimA),
                accum::blitz_existing(elevA)),
                elevAi, false);
        }
    }

    // Get master list of indices we use (in sorted order)
    consolidate(indicesA);
    consolidate(indicesEg);


    // ------ Adjust land type fractions
    for (long iA : indicesA) {
        fgice(iA) = wA(iA) / areaA(iA);
        double fgrnd_val = 1. - (flake(iA) + focean(iA) + fgice(iA));
        fgrnd(iA) = (fgrnd_val < 0 ? 0 : fgrnd_val);

        double focean_val = 1. - (flake(iA) + fgrnd(iA) + fgice(iA));
        if (focean_val < 0) (*icebin_error)(-1,
            "FOCEAN went negative on iA=%ld; take from some other land surface type.", iA);

    }


    // 4) Fix bottom of atmosphere
    for (long iA : indicesA) {
        zatmo(iA) = elevA(iA) * BYGRAV;
    }



    // 3a) Derive FHC from AvE

    // -------- Segment 0: Legacy
    for (long iA : indicesA) {
        fhc2(legacy_segment.base, iA) = 1.e-30;
        underice2(legacy_segment.base, iA) = UI_NOTHING;
        elevE2(legacy_segment.base, iA) = elevA(iA) * fgice(iA);
    }

    // -------- Segment 1: SeaLand
    for (long iA : indicesA) {
        // Sea EC
        fhc2(sealand_segment.base, iA) = 0;
        underice2(sealand_segment.base, iA) = UI_NOTHING;
        elevE2(sealand_segment.base, iA) = 0;

        // Land EC
        fhc2(sealand_segment.base+1, iA) = 1.e-30;
        underice2(sealand_segment.base+1, iA) = UI_NOTHING;
        elevE2(sealand_segment.base+1, iA) = elevA(iA);
    }

    // ------- Segment 2: IceBin
    for (auto ii=begin(*AvE.M); ii != end(*AvE.M); ++ii) {
        iA = dimA.to_sparse(ii->row());
        iE = dimE.to_sparse(ii->col());
        iEg = iE + E_to_Eg;
        auto const iE_tuple(gcm_regridder->indexingHC.index_to_tuple(iE));
            int const iA2 = iE_tuple[0];
            int const ihc = iE_tuple[1];
        if (iA2 != iA) (*icebin_error)(-1,
            "Matrix is non-local: iA=%ld, iE=%ld, iA2=%d", iA, iE, iA2);

        fhcEg(iEg) = ii->value();
        undericeEg(iEg) = UI_ICEBIN;

    }






shapeA = zatmo_m.shape
elevA = elevA1.reshape(shapeA)
areaA = areaA1.reshape(shapeA)
wA = wAvE.reshape(shapeA)
_maskA = ~np.isnan(elevA)
if args.elev_mask_type == 'mar':
    elallA = elallA1.reshape(shapeA)




# ------ Adjust land type fractions
fgice0 = copy.copy(fgice)
fgice[_maskA] = wA[_maskA] / areaA[_maskA]
fgrnd = 1. - (flake + focean + fgice)
fgrnd[fgrnd<0] = 0.
focean = 1. - (flake + fgrnd + fgice)
if not np.all(focean >= 0):
    print('focean', focean[focean<0])
    raise ValueError('FOCEAN went negative; take from some other land surface type')

}





void add_fhc(
    GCMRegridder *gcm_regridder,
    std::vector<std::unique_ptr<DenseVector<1>>> const &elevIs)
{
    for (size_t i=0; i < ice_regridders.size(); ++i) {
        DenseVector<1> const &elevI(*elevIs[i]);
        IceRegridder const &ice_regridder(gcm_regridder->ice_regridders[i]);


        ice_regridder->set_elevI(elevI);
        RegridMatrices rm(ice_regridder);
        

    }
}
