



// ----------------------------------------------------------------
// ----------------------------------------------------------------
/** Diagonal matrix converts wOm <- wOp, weight of the elevation classes in
the two different systems (ModelE (m) vs. IceBin (p) )
NOTES:
 1. focean_m is fixed over the course of a run, since ModelE can not change its ocean mask.
 2. Matrix is diagonal in sparse indexing, not necessarily in dense indexing.
*/
static void raw_OmvOp(
    MageDenseEigenT::AcccumT &&ret,        // {dimA, dimO}
    blitz::Array<double,1> const *_focean_m,    // sparse indexing
    blitz::Array<double,1> const *_focean_p,    // sparse indexing
    blitz::Array<double,1> const *_wOp)        // 0-based, dense indexing
{
    // Make sure we have zero-based indexing
    blitz::Array<double,1> focean_m(reshape1(*_focean_m,0));
    blitz::Array<double,1> focean_p(reshape1(*_focean_p,0));
    blitz::Array<double,1> &wOp(*_wOp);

    SparseIndexT &dimOm(*dims[0]);
    SparseIndexT &dimOp(*dims[1]);

    // By looping over only things in IceBin ice sheet, we implicitly only
    // look at ice on ocean grid cells where both ModelE and IceBin have something.

    // Look only at ocean grid cells where BOTH ModelE and IcBin have ice
    for (int iOp_d=0; iOp_s < nOp_d; ++iOp_d) {
        long const iO_s = dimOp.to_sparse(iOp_d);

        double const fcont_p = 1.0 - focean_p(iO_s);
        double const fcont_m = 1.0 - focean_m(iO_s);

        if (fcont_m == 0.0) continue;
        if (fcont_m != 1.0) (*icebin_error)(-1,
            "fcont_m must be 0 or 1");

        if (fcont_p == 0.0) continue;    // Can't extend ice that's not there

        // In Om space, cell only accounts for the fraction covered by continent.
        // When divided by wM (scaling), this will have the effect of magnifying
        // cells only partially covered by ocean (in PISM)
        ret.add({iO_s, iO_s}, 1./fcont_p);
    }
}
// ----------------------------------------------------------------
static void raw_AvO(
    MageDenseEigenT::AcccumT &&ret,        // {dimA, dimO}
    HntrGrid const &hntrA,
    HntrGrid const &hntrO)
{
    Hntr hntrAO(hntrO, hntrA, 0);    // dimB=A,  dimA=O
    hntrAO.matrix(std::move(ret), 
        std::bind(&dim_clip, &dimA, _1));
}
// ----------------------------------------------------------------
class RawEAvEO {
    MakeDenseEigenT::AccumT ret;    // Final place for EAvEO
    blitz::Array<double,1> const &wEO;

    RawEAvEO(
        MakeDenseEigenT::AccumT &&_ret,
        blitz::Array<double,1> *_wEO)
    : ret(std::move(_ret)), wEO(*_wEO) {}

    /** Called by Hntr::matrix() (AvO) */
    void add(std::array<long,2> index, double value)
    {
        SparseSetT &dimEO(*ret.dims[1]);
        long lA_s = index[0];
        long lO_s = index[1];

        // Iterate through all possible elevation classes for this gridcell pair
        for (int ihc=0; ihc<nhc; ++ihc) {
            long const lEO_s = indexingHCO.tuple_to_index({lO_s,ihc});
            if (!dimEO.in_dense(lEO_s)) continue;    // wEO==0 here
            int const lEO_d = dimEO.to_dense(lEO_s);

            double const weightEO = wEO(lEO_d);
            if (weightEO != 0) {
                int const lEA_s = indexingHCA.tuple_to_index({lA_s,ihc});
                ret.add({lEA_s,lEO_s}, weightEO);
            }
        }
    }
}


static bool clip_true(long ix) { return true; }
/** Can compute EAm_v_EOm or EAi_v_EOi (ModelE or IceBin matrices).
    EA = Elevation grid for A (Atmosphere)
    EO = Elevation grid for O (Ocean)
    Xi = Grid X for IceBin (PISM) ice sheet
    Xm = Grid X for ModelE ice sheet
*/
static void raw_EAvEO(
    MageDenseEigenT::AcccumT &&ret,        // {dimEA, dimEO}
    HntrGrid const *hntrA,
    HntrGrid const *hntrO,
    SparseSetT const *dimA,        // By converting dimO from EOvI
    blitz::Array<double,1> *_wEO)            // == EOvI.wM
{
    Hntr hntrAO(*hntrO, *hntrA, 0);    // dimB=A,  dimA=O
    hntrAO.matrix(
        RawEAvEO(std::move(ret), _wEO),
        std::bind(&dim_clip, &dimA, _1));
}
// ----------------------------------------------------------------
/** Helper function: clip a cell according to its index */
static bool dim_clip(SparseSetT const *dim, long index)
    { return dim->in_sparse(index); }

// ----------------------------------------------------------------
// ========================================================================

/* std::unique_ptr<WeightedSparse> compute_EAmvIp(
    std::array<SparseSet *,2> dims,
    RegridMatrices::Params const &paramsA,
    IceRegridder_ModelE *regridderA,
    RegridMatrices const &rmO)
*/
#define URAE E
#include "compute_EA.hpp"
#undef URAE

#define URAE A
#include "compute_EA.hpp"
#undef URAE
// ========================================================================

// ========================================================================
/** Creates Atmosphere grid from an existing Hntr-type Ocean grid. */
static std::unique_ptr<Grid> make_gridA(Grid const *_gridO)
{
    // -------- Check types on gridO
    Grid_LonLat const *gridO(dynamc_cast<Grid_LonLat *>(_gridO));
    if (!gridO) (*icebin_error)(-1,
        "make_gridA() requires type Grid_LonLat");

    if (gridO->north_pole != gridO->south_pole) (*icebin_error)(-1,
        "north_pole=%d and south_pole=%d must match in gridO",
        gridO->north_pole, gridO->south_pole);

    HntrGrid const *_hntrO(&*gridO->hntr);
    if (!_hntrO) (*icebin_error)(-1,
        "make_gridA() requires gridO have a Hntr source");
    HntrGrid const &hntrO(*_hntrO);

    if ((hntrO.im % 2 != 0) || (hntrO.jm % 2 != 0)) (*icebin_error)(-1,
        "Ocean grid must have even number of gridcells for im and jm (vs. %d %d)",
        hntrO.im, hntrO.jm);

    // --------------------------
    // Define Atmosphere grid to be exactly twice the Ocean grid
    HntrGrid const hntrA(hntrO.im/2, hntrO.jm/2, hntrO.offi*0.5, hntrO.dlat*2.);

    // Use hntr to figure out which grid cells should be realized in A, based
    // on realized grid cells in O
    SparseSetT dimA;
    Hntr hntrOvA(hntrA, hntrO, 0);    // gridB=O, gridA=A
    hntrOvA.matrix(
        accum::SparseSetAccum<SparseSetT,double,2>({nullptr, &dimA}),
        std::bind(&dim_clip, &gridO->dim, _1));

    // ------- Produce a full gridA, based on hntrA and realized cells in dimA
    GridSpec_Hntr spec(hntrA);
    spec.name = gridO->name + "_A";
    spec.pole_caps = gridO->north_pole;

    // Keep cells listed in dimA
    spec.spherical_clip = std::bind(&dim_clip, &dimA, _1);
    spec.points_in_side = 1;    // Keep it simple, we don't need this anyway
    spec.eq_rad = gridO->eq_rad;

    std::unique_ptr<Grid_LonLat *> gridA(new Grid_LonLat);
    spec.make_grid(*gridA);
}


GCMRegridder_ModelE(std::unique_ptr<icebin::GCMRegridder> &&_gcmO) :
    : gcmO(std::move(_gcmO))
{
    // Initialize baseclass members
    gridA = make_gridA(gcmO->gridA);
    correctA = gcmO->correctA;
    hcdefs = gcmO->hcdefs;
    indexingHC = gcmO->indexingHC;
    indexingE = derive_indexingE(gridA->indexing, indexingHC);
}

RegridMatrices const GCMRegridder_ModelE::regrid_matrices(std::string const &ice_sheet_name) const
{
    IceModelRegridder const *ice_regridderO = gcmO->ice_regridder(ice_sheet_name);

    // Navigate the data structures
    Grid_LonLat const *gridO = dynamic_cast<Grid_LonLat const *>(&*gcmO->gridA);
    HntrGrid const &hntrO(*gridO->hntr);

    // Construct set of cells in A
    SparseSetT dimA;
    for (auto cell = gridA->cells.begin(); cell != gridA->cells.end(); ++cell)
        dimA.add_dense(cell->index);

    RegridMatrices rm;
    RegridMatrices const &rmO(rm.tmp.make<RegridMatrices>(
        ice_regridderO->regrid_matrices()));
    rm.add_regrid("EAmvIp", std::bind(&compute_EAmvIp, _1, _2,
        &dimA, &rmO));
    rm.add_regrid("AmvIp", std::bind(&compute_AmvIp, _1, _2,
        &dimA, &rmO));

    return rm;
}
