class ConstUniverse {
    std::vector<std::string> names;
    std::vector<SparseSetT *> dims;
    std::vector<int> extents;

    ConstUniverse(
        std::vector<std::string> &&_names
        std::vector<SparseSetT *> &&_dims) :
        names(std::move(_names)), dims(std::move(_dims))
    {
        if (names.size() != dims.size()) (*icebin_error(-1,
            "names.size() and dims.size() must match"));

        extents.reserve(dims.size());
        for (size_t i=0; i<dims.size(); ++i)
            extents.push_back(dims[i]->dense_extent());
    }

    ~ConstUniverse()
    {
        bool err = false;
        for (size_t i=0; i<dims.size(); ++i) {
            if (extents[i] != dims[i]->dense_extent()) {
                fprintf(stderr, "Dimension %s changed from %d to %d\n",
                    names[i].c_str(), extents[i], dims[i]->dense_extent());
                err = true;
            }
        }
        if (err) (*icebin_error(-1,
            "At least one dimension changed"));
    }
};


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
    SparseSetT const &dimXAm,
    RegridMatrices const &rmO)
{
    TmpAlloc tmp;

    SparseSetT &dimXAm(*dims[0]);
    SparseSetT &dimIp(*dims[1]);
    SparseSetT dimEOm, dimEOp, dimOm;

    std::unique_ptr<WeightedSparse> ret(new WeightedSparse(dims, false));    // not conservative

    // ------------ Params for generating sub-matrices
    RegridMatrices::Params paramsO(paramsA);
    paramsO.scale = false;





    // ---------------- Take along the non-conservative weights



    // ----------- Determine universe for AOp / EOp
    // (and we will use the main matrix later)

    // We need AOpvIp for wAOP; used below.
    SparseSetT dimAOp;
    std::unique_ptr<WeightedSparse> AOpvIp(rmO.matrix("AvI", {&dimAOp, &dimIp}, paramsO));
    blitz::Array<double,1> const &wAOp(AOpvIp.mW);

    // dimAOm is a subset of dimAOp; we will use dimAOp to be sure.
    SparseSetT &dimAOm(dimAOp);

    std::unique_ptr<ConstUniverse> const_universe;

#if STR(URAE) == "E"

    // We also need to get the universe for EOp / EOm
    SparseSetT dimEOp;
    const_universe.reset(new ConstUniverse({"dimIp"}, {&dimIp});
    std::unique_ptr<WeightedSparse> EOpvIp(rmO.matrix("EvI", {&dimEOp, &dimIp}, paramsO));
    const_universe.reset();        // Check that dimIp hasn't changed

    // dimEOm is a subset of dimEOp
    SparseSetT &dimEOm(dimEOp);
#endif

    // ----------- Compute weights for GCM ice sheet (wXOm)
    // wEOm = sEOmvOm * EOmvAOm * wc_AOmvAOp * wAOp
    //    or   wAOm = AOmvAOp * wAOp
    // wAOp

    // OmvOp
    EigenSparseMatrixT sc_AOmvAOp(MakeDenseEigenT(
        std::bind(&weightconv_AOmvAOp, _1, focean_m, focean_p, &wAOp),
        {SparsifyTransform::TO_DENSE},
        {&dimAOm, &dimAOp}, '.').to_eigen());

    EigenColVectorT wAOm_e(wc_AOmvAOp * map_eigen_colvector(wAOp));


    EigenColVectorT *wXOm_e;
    if (X == 'E') {
        // EOmvAOm: Repeat A values in E
        const_universe.reset(new ConstUniverse({"dimEOm", "dimAOm"}, {&dimEOm, &dimAOm});
        WeightedSparse EOmvAOm(rmO.matrix("EvA", {&dimEOm, &dimAOm}, paramsO));
        const_universe.reset();        // Check that dims didn't change
        auto sEOmvOm(1. / EOmvAOm->wM);

        // wEOm_e
        tmp.reset_ptr(wXOm_e, map_eigen_diagonal(sEOmvAOm) * *EOmvOm->M * wAOm_e);
    } else {
        wXOm_e = &wAOm_e;
    }

    // ---------------- Compute the main matrix
    // ---------------- XAmvIp = XAmvXOm * XOpvIp
    std::unique_ptr<EigenSparseMatrixT> XAmvXOm;    // Unscaled matrix
    if (X == 'E') {
        SparseSetT &dimEAm(dimXAm);
        reset_ptr(XAmvXOm, MakeDenseEigenT(
            std::bind(&raw_EAvEO, _1, &hntrA, &hntrO, &dimEOm, &wEOm),
            {SparsifyTransform::ADD_DENSE, SparsifyTransform::TO_DENSE_IGNORE_MISSING},
            {&dimEAm, &dimEOm}, '.').to_eigen());
    } else {
        SparseSetT &dimAAm(dimXAm);
        // Actually AAmvAOm
        reset_ptr(XAmvXOm, MakeDenseEigenT(
            std::bind(&raw_AAvAO, _1, &hntrA, &hntrO),
            {SparsifyTransform::ADD_DENSE},
            {&dimAAm, &dimAOm}, '.').to_eigen());
    }

    // ---------- wEAm = [scale] * XAmvXOm * wEOm
    blitz::Array<double,1> sXAmvXOm(sum(XAmvXOm, 0, '-'));
    auto &wXAm_e(ret->tmp.make<EigenColVectorT>(
        map_eigen_diagonal(sXAmvXOm) * *XAmvXOm * *wXOm_e));

    // ----------- Put it all together (XAmvIp)
    ret->wM.reference(to_blitz(wXAm_e));
    if (paramsA.scale) {
        blitz::Array<double,1> sXAm(1. / to_blitz(wXAm_e));
        ret->M.reset(new EigenSparseMatrixT(
            map_eigen_diagonal(sXAm) * *XAmvXOm * *XOpvIp->M));
    } else {
        ret->M.reset(new EigenSparseMatrixT(
            *XAmvXOm * *XOpvIp->M));
    }
    sum->Mw.reference(EOpvIp->Mw);
    return ret;
}

#undef STR
