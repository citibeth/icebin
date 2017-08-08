

struct DualiceSparse : public WeightedSparse {
    blitz::Array<double,1> &wM_p;
    blitz::Array<double,1> &Mw_p;
    blitz::Array<double,1> wM_m;
    blitz::Array<double,1> Mw_m;

    DualiceSparse(WeightedSparse const &base) :
        WeightedSparse(base.dims),
        wM_p(wM), Mw_p(Mw),
        wM_m(base.dims[0]->dense_extent()),
        Mw_m(base.dims[0]->dense_extent())
    {
        wM_p.reference(blitz::Array<double,1>(OvE_p.dims[0]->dense_extent());
        wM_m.reference(blitz::Array<double,1>(OvE_p.dims[0]->dense_extent());

        Mw_p.reference(blitz::Array<double,1>(OvE_p.dims[1]->dense_extent());
        Mw_m.reference(blitz::Array<double,1>(OvE_p.dims[1]->dense_extent());

        smooth = OvE_p.smooth;
        conserve = OvE_p.conserve;
    }
};

/** Makes an OvE matrix for ModelE ice sheet, based on the same for the PISM ice sheet.
Before this is called, make sure that dimO and dimE are initialized with all places where focean_p != 0 or focean_m != 0.  Works on scaled or unscaled OvE_p. */
std::unique_ptr<WeightedSparse> new_OvE_m(
    WeightedSparse const &OvE_p,
    blitz::Array<double,1> const &focean_p,    // 0-based
    blitz::Array<double,1> const &focean_m)    // 0-based
{
    if ((focean_p.lbound(0) != 0) || (focean_m.lbound(0) != 0)) (*icebin_error)(-1,
        "focean arrays must be zero-based");

    int const _dimO = 0;
    int const _dim1 = 1;
    auto &dimO(OvE_p.dims[0]);
    auto &dimE(OvE_p.dims[1]);
    std::unique_ptr<DualiceSparse> OvE_m(OvE_p);

    // Compute factor by which basis functions enlarge in cells of M vs. P
    int const nO = dimO.dense_extent();
    int const nE = dimE.dense_extent();
    for (int iO=0; iO<nO; ++iO) {
        double const fcont_p = 1.0 - focean_p(iO);
        double const fcont_m = 1.0 - focean_m(iO);

        if (fcont_p == fcont_m) {
            OvE_m->wM(iO) = OvE_p.wM(iO);
        } else if (fcont_m == 1.) {    // (&& 0 < fcont_p < 1)
            // Everything in this grid cell gets expanded in ModelE ice sheet
            OvE_m->wM(iO) = OvE_p.wM(iO) / fcont_p;
        } else {
            // This cell not used by ModelE ice sheet
            OvE_m->wM(iO) = 0;
        }
    }

    // 1. Convert mW (in O space) to Wm (in E space)
    // 2. Convert main matrix
    blitz::Array<double,1> mbypE(dimE.dense_extent());
    TupleList<int,double,2> M;
    for (auto ii(begin(*OvE_p.M)); ii != end(*OvE_p.M); ++ii) {
        auto const iO(ii->row());
        auto const iE(ii->col());
        E2O(iE) = iO;

        double const fcont_p = 1.0 - focean_p(iO);
        double const fcont_m = 1.0 - focean_m(iO);

        OvE_m->Mw(iE) = OvE_m.wM(iO);
        double const mbyp = OvE_m.wM(iO) / OvE_p.wM(iO);
        if (mbyp != 0) {
            M.add({iO,iE}, ii->value() * mbyp);
        }
    }
    OvE_m->M.reset(new EigenSparseMatrixT(to_eigen(M)));
    return OvE_m;
}

static void nop(MakeDenseEigenT::AccumT &accum) {}

std::unique_ptr<WeightedSparse> AvX_v_OvX(
    WeightedSparse const &OvX,
    Indexing const &indexingO,
    Hntr const &hntr_AvO)    // Converts A <- O
{
    std::unique_ptr<WeightedSparse> AvX(new WeightedSparse(OvX));

    SparseSet const &dimO(*OvX->dims[0]);
    SparseSet const &dimX(*OvX->dims[1]);

    // Initialize dimA from dimO
    SparseSet &dimA(AvX.tmp.make<SparseSet>());
    for (int iO_d=0; iO_d<dimO.dense_extent(); ++iO_d) {
        long iO_s = dimO.to_sparse(iO_d);

        std::array<int,2> tO(indexingO.index_to_tuple(iO_s));
        std::array<int,2> tA0 { tO[0]*2, tO[1]*2 };
        for (delta_j=0; delta_j<2; ++delta_j) {
        for (delta_i=0; delta_i<2; ++delta_i) {
            int iA_s = tuple_to_index(std::array<int,2>
                { tA0[0] + delta_i, tA0[1] + delta_j });
            dimA.add(iA_s);
        }}
    }

    // Get the sparse matrix to convert A<-O in sparse (complete) indexing
    spsparse::TupleList<int,double,2> &AvO;
    HntrGrid const &Ogrid(hntr_AvO.Agrid);
    blitz::Array<double, 1> WTO(const_array(    // 1-based
        blitz::shape(Ogrd.size()), 1.0, FortranArray<1>()));
    hntr_AvO.matrix(AvO, WTO);    // AvO is 0-based; WTO 1-based

    // Get the AvO matrix (just the subset we need)
    EigenSparseMatrixT AvO(MakeDenseEigenT(
        std::bind(&Hntr::matrix, std::placeholders::_1,
            const_array(blitz::shape(Ogrid.size()), 1.0, FortranArray<1>())),
        SparsifyTransform::TO_DENSE_IGNORE_MISSING,
        {&dimA, &dimX}, '.').to_eigen());

    // Use the AvO matrix to convert AvX
    AvX->wM = AvO * OvX->wM;
    AvX->M.reset(new EigenSparseMatrixT(AvO * *AvX->M));
    AvX->Mw = OvX->Mw;

    return AvX;
}


/** Computes I_p <-- E_m
Eg: converts PISM ice sheet (Ice grid) to ModelE ice sheet (Elevation grid) */
EigenDenseMatrixT apply_dualice_IvE_e(
    WeightedSparse const &IvE_p,
    blitz::Array<double,2> const &E_b,    // nvar x nE: The fields we want to convert
    blitz::Array<double,1> const &wE_p_b,    // From AvE_p
    blitz::Array<double,1> const &wE_m_b,    // From AvE_m
    double fill)
{
    // Mass of E_p and E_m basis vectors   [1 x nE]
    Eigen::Map<EigenRowVectorT> const wE_p(
        const_cast<double *>(wE_p_b.data()),
        wE_p_b.extent(0));
    Eigen::Map<EigenRowVectorT> const wE_m(
        const_cast<double *>(wE_m_b.data()),
        wE_m_b.extent(0));

    // E_b --> Eigen Matrix   [nE x nvar]
    Eigen::Map<EigenDenseMatrixT> const E(
        const_cast<double *>(E_b.data()),
        E_b.extent(1), E_b.extent(0));

    // Correct each variable by this ratio [nvar x nvar]
    auto correct( ((wE_m*E).array() / (wE_p*E).array()) .asDiagonal());
    return IvE_p.apply_e(E * correct, fill); * ;
}


EigenDenseMatrixT apply_dualice_AvI_e(
    WeightedSparse const &AvI_p,
    blitz::Array<double,2> const &I_b,    // nvar x nI: The fields we want to convert
    blitz::Array<double,1> const &wA_p_b,    // From AvE_p
    blitz::Array<double,1> const &wA_m_b,    // From AvE_m
    double fill)
{
    // Mass of E_p and E_m basis vectors   [1 x nE]
    Eigen::Map<EigenRowVectorT> const wA_p(
        const_cast<double *>(wA_p_b.data()),
        wA_p_b.extent(0));
    Eigen::Map<EigenRowVectorT> const wA_m(
        const_cast<double *>(wA_m_b.data()),
        wA_m_b.extent(0));

    // I_b --> Eigen Matrix   [nE x nvar]
    Eigen::Map<EigenDenseMatrixT> const I(
        const_cast<double *>(I_b.data()),
        I_b.extent(1), I_b.extent(0));

    // Correct each variable by this ratio [nvar x nvar]
    auto correct( ((wA_p*E).array() / (wA_m*E).array()) .asDiagonal());
    return AvI_p.apply_e(I, fill) * correct;
}


