// Included as part of IceRegridder_ModelE.cpp
// Exploits symmetry in matrix construction between A and E grids.
// Using the C++ preprocessor to make two versions was deemed more
// clear in this case than using templates or std::function or
// something else.
//
// The symbol URAE must be declared to be either:

#define STR(x) #x

std::unique_ptr<WeightedSparse> compute_ ## URAE ## mvIp(
    std::array<SparseSet *,2> dims,
    RegridMatrices::Params const &paramsA,
    SparseSetT const &dimA,
    RegridMatrices const &rmO)
{
    SparseSetT &dimEAm(dims[0]);
    SparseSetT &dimIp(dims[1]);
    SparseSetT dimEOm, dimEOp, dimIp, dimOm;

    // ------------ Params for generating sub-matrices
    RegridMatrices::Params paramsO(paramsA);
    paramsO.scale = false;


    // ---------------- Take along the non-conservative weights

    // ----------- wEOm = sEOmvOm * EOmvOm * OmvOp * wOp
    // wOp
    auto OpvIp(rmO.matrix("AvI"));
    blitz::Array<double,1> const &wOp(OpvIp.mW);

    // OmvOp
    auto OmvOp(MakeDenseEigenT(
        std::bind(&raw_OmvOp, _1, focean_m, focean_p, &wOp),
        {SparsifyTransform::ADD_DENSE},
        {&dimOm, &dimOp}, '.').to_eigen());

#if STR(URAE) == "E"
    // EOmvOm
    auto EOmvOm(rmO.matrix("EvA", {&dimEOm, &dimOm}, paramsO));
    auto sEOmvOm(1. / EOmvOm->wM);

    EigenColVectorT wEOm_e(
        map_eigen_diagonal(sEOmvOm) * EOmvOm->M * OmvOp * map_eigen_colvector(wOp));
    blitz::Array<double,1> wEOm(to_blitz(wEOm_e));
#else
    // Actually wOm
    blitz::Array<double,1> wEOm(OmvOp->wM);
#endif

    // ---------------- Compute the main matrix
    // ---------------- EAmvIp = EAmvEOm * EOpvIp

    std::unique_ptr<WeightedSparse> EOpvIp(rmO.matrix("EvI"));

#if STR(URAE) == "E"
    auto EAmvEOm(MakeDenseEigenT(
        std::bind(&raw_EAvEO, _1, &hntrA, &hntrO, &dimA, &wEOm),
        {SparsifyTransform::ADD_DENSE},
        {&dimEAm, &dimEOm}, '.').to_eigen());
#else
    // Actually AmvOm
    auto EAmvEOm(MakeDenseEigenT(
        std::bind(&raw_AvO, _1, &hntrA, &hntrO),
        {SparsifyTransform::ADD_DENSE},
        {&dimEAm, &dimEOm}, '.').to_eigen());
#endif

    // ---------- wEAm = EAmvEOm * wEOm
    auto wEAm(EAmvEOm * wEOm_e);

    // ----------- Put it all together (EAmvIp)
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse(dims, true));
    ret->wM.reference(wEAm);
    if (paramsA.scale) {
        blitz::Array<double,1> sEAm(1. / wEAm);
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(sEAm) * EAmvEOm * *EOpvIp->M));
    } else {
        ret->M.reset(new EigenSparseMatrixT(
            EAmvEOm * EOpvIp));
    }
    sum->Mw.reference(EOpvIp.Mw);
    return ret;
}

#undef STR
