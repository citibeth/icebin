/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include <unordered_set>

#include <spsparse/netcdf.hpp>
#include <spsparse/multiply_sparse.hpp>

#include <icebin/GCMRegridder.hpp>
#include <icebin/IceRegridder_L0.hpp>

using namespace std;
using namespace netCDF;
using namespace ibmisc;
using namespace std::placeholders;  // for _1, _2, _3...
using namespace spsparse;

namespace icebin {

// -----------------------------------------------------
IceRegridder::IceRegridder() : interp_style(InterpStyle::Z_INTERP), _name("icesheet") {}

IceRegridder::~IceRegridder() {}


std::unordered_map<long,double> IceRegridder::elevI_hash() const
{
    // Convert elevI to a hash table, so we can look up in it easily
    std::unordered_map<long,double> elevIh;
    for (auto ii=elevI.begin(); ii != elevI.end(); ++ii)
        elevIh.insert(std::make_pair(ii.index(0), ii.val()));
    return elevIh;
}
// --------------------------------------------------------

// ==============================================================
// Different vector spaces:
//      Description                    Domain
// ---------------------------------------------
// A  = Atmosphere grid                sphere
// Ap = Projected atmosphere grid      plane
// E  = Elevation grid                 sphere
// Ep = Projected elevation grid       plane
// I  = Ice grid                       plane
//
// Write out the parts that this class computed --- so we can test/check them
// -------------------------------------------------------------
/** Produces the diagonal matrix [Atmosphere projected] <-- [Atmosphere]
NOTE: wAvAp == sApvA */
void IceRegridder::sApvA(SparseTriplets<SparseMatrix> &w, std::function<bool(long)> const &filter_fn)
{
    ibmisc::Proj_LL2XY proj(gridI->sproj);
    for (auto cell=gcm->gridA->cells.begin(); cell != gcm->gridA->cells.end(); ++cell) {
        if (!filter_fn(cell->index)) continue;

        w.add({cell->index, cell->index}, cell->native_area / cell->proj_area(&proj));
    }
}

/** Produces the diagonal matrix [Elevation projected] <-- [Elevation]
NOTE: wAvAp == sApvA */
void IceRegridder::sEpvE(SparseTriplets<SparseMatrix> &w, std::function<bool(long)> const &filter_fn)
{
    ibmisc::Proj_LL2XY proj(gridI->sproj);
    for (auto cell=gcm->gridA->cells.begin(); cell != gcm->gridA->cells.end(); ++cell) {
        long nhp = gcm->nhp(cell->index);
        long tuple[2] = {cell->index, 0};
        long &ihp(tuple[1]);
        for (ihp=0; ihp<nhp; ++ihp) {
            long indexE = gcm->indexingHC.tuple_to_index(tuple);
            if (!filter_fn(indexE)) continue;
            w.add({indexE, indexE}, cell->native_area / cell->proj_area(&proj));
        }
    }
}
// -------------------------------------------------------------
void IceRegridder::clear()
{
    gridI.reset();
    exgrid.reset();
    elevI.clear();
}
// -------------------------------------------------------------
/** Determines correct indexing for the E grid.
Returns: im,jm,...,iHC*/
static ibmisc::Indexing<int,long> derive_indexingE(
ibmisc::Indexing<int,long> const &indexingA,      // im,jm,...
ibmisc::Indexing<long,long> const &indexingHC)    // iA,iHC
{
    std::vector<std::string> name(indexingA.name);    // Name of each dimension
    std::vector<int> base(indexingA.base);   // First element in each index
    std::vector<int> extent(indexingA.extent); // Extent (# elements) of each index
    std::vector<int> indices(indexingA.indices);   // Index IDs sorted by descending stride. {0,1,...} for row-major, reversed for col-major

    name.push_back(indexingHC.name[1]);
    base.push_back(indexingHC.base[1]);
    extent.push_back(indexingHC.extent[1]);

    if (indexingHC.indices[0] == 0) {    // HC is row major so we are too
        indices.push_back(indexingA.rank());
    } else {
        for (auto &ii : indices) ++ii;
        indices.push_back(0);
    }

    return ibmisc::Indexing<int,long>(
        std::move(name), std::move(base),
        std::move(extent), std::move(indices));
}
// -------------------------------------------------------------
void IceRegridder::ncio(NcIO &ncio, std::string const &vname)
{
    if (ncio.rw == 'r') {
        clear();
        gridI = new_grid(ncio, vname + ".gridI");
        exgrid = new_grid(ncio, vname + ".exgrid");
    }

    auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});
    get_or_put_att(info_v, ncio.rw, "name", _name);
    get_or_put_att_enum(info_v, ncio.rw, "interp_style", interp_style);

    gridI->ncio(ncio, vname + ".gridI");
    exgrid->ncio(ncio, vname + ".exgrid");
    ncio_spsparse(ncio, elevI, true, vname + ".elevI");

    indexingE = derive_indexingE(gridA->indexing, indexingHC);
}
// ========================================================
void GCMRegridder::init(
    std::unique_ptr<Grid> &&_gridA,
//  ibmisc::Domain<int> &&_domainA,     // Tells us which cells in gridA to keep...
    std::vector<double> &&_hcdefs,  // [nhp]
    ibmisc::Indexing<long,long> &&_indexingHC,
    bool _correctA)
{
    gridA = std::move(_gridA);
//  domainA = std::move(_domainA);
    hcdefs = std::move(_hcdefs);
    indexingHC = std::move(_indexingHC);
    correctA = _correctA;

    if (indexingHC.rank() != 2) (*icebin_error)(-1,
        "indexingHC has rank %d, it must have rank=2", indexingHC.rank());

    indexingE = derive_indexingE(gridA->indexing, indexingHC);
}
// -------------------------------------------------------------
std::unique_ptr<IceRegridder> new_ice_regridder(IceRegridder::Type type)
{
    switch(type.index()) {
        case IceRegridder::Type::L0 :
            return std::unique_ptr<IceRegridder>(new IceRegridder_L0);
        break;
        default :
            (*icebin_error)(-1,
                "Unknown IceRegridder::Type %s", type.str());
        break;
    }
}

std::unique_ptr<IceRegridder> new_ice_regridder(NcIO &ncio, std::string vname)
{
    std::string vn(vname + ".gridI.info");
    auto gridI_info_v = get_or_add_var(ncio, vn, "int64", {});

    Grid::Parameterization parameterization;
    get_or_put_att_enum(gridI_info_v, ncio.rw, "parameterization", parameterization);
    return new_ice_regridder(parameterization);
}

// -------------------------------------------------------------
// ==============================================================

void GCMRegridder::clear()
{
    gridA.reset();
    hcdefs.clear();
    sheets_index.clear();
    sheets.clear();
}

void GCMRegridder::ncio(NcIO &ncio, std::string const &vname)
{
    auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});

    if (ncio.rw == 'r') {
        clear();
        gridA = new_grid(ncio, vname + ".gridA");   // Instantiates but does not read/write
    }

    // Read/Write gridA and other global stuff
    gridA->ncio(ncio, vname + ".gridA");
    indexingHC.ncio(ncio, ncInt, vname + ".indexingHC");
    ncio_vector(ncio, hcdefs, true, vname + ".hcdefs", ncDouble,
        get_or_add_dims(ncio, {vname + ".nhp"}, {hcdefs.size()} ));
//  domainA.ncio(ncio, ncInt, vname + ".domainA");
    get_or_put_att(info_v, ncio.rw, vname + ".correctA", &correctA, 1);

    // Read/write list of sheet names
    std::vector<std::string> sheet_names;
    if (ncio.rw == 'w') {
        // Copy sheet names to a std::vector
        for (auto ii=sheets_index.begin(); ii != sheets_index.end(); ++ii)
            sheet_names.push_back(*ii);
    }
    get_or_put_att(info_v, ncio.rw, vname + ".sheets", sheet_names);

    // Read/write the regridder for each ice sheet
    if (ncio.rw == 'r') {   // Instantiate
        for (auto sheet_name = sheet_names.begin(); sheet_name != sheet_names.end(); ++sheet_name) {
            std::string vn(vname + "." + *sheet_name);
            add_sheet(*sheet_name, new_ice_regridder(ncio, vn));
        }
    }
    for (auto sheet=sheets.begin(); sheet != sheets.end(); ++sheet) {
        (*sheet)->ncio(ncio, vname + "." + (*sheet)->name());
    }
}
// -------------------------------------------------------------

// ---------------------------------------------------------------------
/** Made for binding... */
static bool in_good(std::unordered_set<long> const *set, long index_c)
{
    return (set->find(index_c) != set->end());
}

void IceRegridder::filter_cellsA(std::function<bool (long)> const &useA)
{

  // Figure out which cells to keep

    // List of cells in gridI / exgrid that overlap a cell we want to keep
    std::unordered_set<long> good_index_gridI;
    std::unordered_set<long> good_index_exgrid;


    std::unordered_set<int> good_j;
    for (auto excell = exgrid->cells.begin(); excell != exgrid->cells.end(); ++excell) {
        int index1 = excell->i;
        if (useA(index1)) {
            good_index_gridI.insert(excell->j);
            good_index_exgrid.insert(excell->index);
        }
    }

    // Remove unneeded cells from gridI
    gridI->filter_cells(std::bind(&in_good, &good_index_gridI, _1));
    exgrid->filter_cells(std::bind(&in_good, &good_index_exgrid, _1));
}

void GCMRegridder::filter_cellsA(std::function<bool(long)> const &keepA)
{

    // Now remove cells from the exgrids and gridIs that
    // do not interact with the cells we've kept in grid1.
    for (auto sheet=sheets.begin(); sheet != sheets.end(); ++sheet) {
        (*sheet)->filter_cellsA(keepA);
    }

    gridA->filter_cells(keepA);
}

void GCMRegridder::filter_cellsA(ibmisc::Domain<int> const &domainA)
{
    filter_cellsA(std::bind(&ibmisc::in_domain<int, long>,
        &domainA, &gridA->indexing, _1));
}

void GCMRegridder::wA(SparseVector &w) const
{
    for (auto cell=gridA->cells.begin(); cell != gridA->cells.end(); ++cell)
        w.add({cell->index}, cell->native_area);
}
// ---------------------------------------------------------------------
void IceRegridder::init(
    std::string const &name,
    std::unique_ptr<Grid> &&_gridI,
    std::unique_ptr<Grid> &&_exgrid,
    InterpStyle _interp_style,
    SparseVector &&_elevI)
{
    _name = (name != "" ? name : gridI->name);
    gridI = std::move(_gridI);
    exgrid = std::move(_exgrid);
    interp_style = _interp_style;
    elevI = std::move(_elevI);
}
// ==============================================================
extern void linterp_1d(
    std::vector<double> const &xpoints,
    double xx,
    int *indices, double *weights)  // Size-2 arrays
{
    int n = xpoints.size();


    // This is the point ABOVE our value.
    // (i0 = i1 - 1, xpoints[i0] < xx <= xpoints[i1])
    // See: http://www.cplusplus.com/reference/algorithm/lower_bound/
    int i1 = lower_bound(xpoints.begin(), xpoints.end(), xx) - xpoints.begin();

    if (i1 <= 0) i1 = 1;
    if (i1 >= n) i1 = n-1;

    int i0 = i1-1;
    indices[0] = i0;
    indices[1] = i1;
    double ratio = (xx - xpoints[i0]) / (xpoints[i1] - xpoints[i0]);
    weights[0] = (1.0 - ratio);
    weights[1] = ratio;
}
// ========================================================================

/** Used to generate Ur matrices for Atm or Elevation grids.
This allows us to exploit an algebraic symmetry between A and E. */
struct UrAE {
    const long nfull;

    typedef std::function<void(spsparse::SparseTriplets<SparseMatrix> &)> GvAp_fn;
    const GvAp_fn GvAp;

    typedef const std::function<void(
        spsparse::SparseTriplets<SparseMatrix> &w,
        std::function<bool(long)>
    )> sApvA_fn;
    const sApvA_fn sApvA;

    UrAE(long _nfull, GvAp_fn _GvAp, sApvA_fn _sApvA) : nfull(_nfull), GvAp(_GvAp), sApvA(_sApvA) {}
};



// ------------------------------------------------------------
std::function<bool(long)> in_sparse_fn(SparseSetT const &dim)
    { return std::bind(&SparseSetT::in_sparse, &dim, _1); }
// ------------------------------------------------------------
static std::unique_ptr<WeightedSparse> compute_AEvI(IceRegridder *sheet, bool scale, bool correctA, UrAE const &AE)
{
    SparseSetT dimA, dimG, dimI;

    // ----- Get the Ur matrices (which determines our dense dimensions)
    // GvI
    SparseTriplets<SparseMatrix> GvI_t({&dimG, &dimI});    // _t=triplets
    GvI_t.set_shape({sheet->nG(), sheet->nI()});
    sheet->GvI(GvI_t);

    // ApvG
    SparseTriplets<SparseMatrix> GvAp_t({&dimG, &dimA});
    GvAp_t.set_shape({sheet->nG(), AE.nfull});
    AE.GvAp(GvAp_t);

    // ----- Convert to Eigen and multiply
    auto GvI_e(GvI_t.to_eigen());
    auto sGvI_e(scale_matrix(sum(GvI_e,0)))
    auto ApvG_e(GvAp_t.to_eigen('T'));
    std::unique_ptr<EigenSparseMatrix> ApvI_e(
        new EigenSparseMatrix(ApvG_e * sGvI_e * GvI_e));

    // ----- Apply final scaling, and convert back to sparse dimension
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse);

    if (correctA) {

        // Scaling matrix AvAp
        SparseTriplets<SparseMatrix> sApvA({&dimA, &dimA});
        sApvA.set_shape({AE.nfull, AE.nfull});
        AE.sApvA(sApvA, in_sparse_fn(dimA));   // Filter to avoid adding more items to dimA

        if (scale) {
            // Get two diagonal Eigen scale matrices
            auto sAvAp_e(sApvA.to_eigen('.', true));        // Invert it...
            auto sApvI_e(weight_matrix(*ApvI_e, 0, true));

            ret->M.reset(new EigenSparseMatrix(
                sAvAp_e * sApvI_e * *ApvI_e));    // AvI_scaled
        } else {
            ret->M = std::move(ApvI_e);
        }

        // ----- Compute the final weight matrix
        auto wAvAp_e(sApvA.to_eigen('.', false));       // Invert it...
        auto wApvI_e(weight_matrix(ApvI_e, 0, false));
        EigenSparseMatrix weight_e(wAvAp_e * wApvI_e);

        // Sums rows of an Eigen matrix into a dense blitz::Array
        ret->weight.reference(sum_rows(weight_e));
    } else {

        if (scale) {
            // Get two diagonal Eigen scale matrices
            auto sApvI_e(weight_matrix(ApvI_e, 0, true));

            ret->M.reset(new EigenSparseMatrix(
                sApvI_e * ApvI_e));    // ApvI_scaled
        } else {
            ret->M = std::move(ApvI_e);
        }

        // ----- Compute the final weight matrix
        auto wApvI_e(weight_matrix(ApvI_e, 0, false));
        ret->weight.reference(sum_rows(wApvI_e));
    }
    ret->dim[0] = std::move(dimA);
    ret->dim[1] = std::move(dimI);

    return ret;
}

static std::unique_ptr<WeightedSparse> compute_IvAE(IceRegridder *sheet, bool scale, bool correctA, UrAE const &AE)
{
    SparseSetT dimA, dimG, dimI;

    // ----- Get the Ur matrices (which determines our dense dimensions)
    // ApvG
    SparseTriplets<SparseMatrix> GvAp({&dimG, &dimA});
    GvAp.set_shape({sheet->nG(), AE.nfull});
    AE.GvAp(GvAp);

    // GvI
    SparseTriplets<SparseMatrix> GvI({&dimG, &dimI});
    GvI.set_shape({sheet->nG(), sheet->nI()});
    sheet->GvI(GvI);

    // ----- Convert to Eigen and multiply
    auto GvAp_e(GvAp.to_eigen('.'));
    auto sGvAp_e(GvAp.eigen_scale_matrix(0));
    auto IvG_e(GvI.to_eigen('T'));

    // Unscaled matrix
    EigenSparseMatrix IvAp_e(IvG_e * sGvAp_e * GvAp_e);

    // ----- Apply final scaling, and convert back to sparse dimension
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse({dimI.sparse_extent(), dimA.sparse_extent()}));
    if (correctA) {
        // sApvA
        SparseTriplets<SparseMatrix> sApvA({&dimA, &dimA});
        sApvA.set_shape({AE.nfull, AE.nfull});
        AE.sApvA(sApvA, in_sparse_fn(dimA)); // Filter to avoid adding more items to dimA
        auto sApvA_e(sApvA.to_eigen('.', false));       // Don't invert it...

        if (scale) {
            auto sIvAp_e(weight_matrix(IvAp_e, 0, true));
            ret->M = sIvAp_e * IvAp_e * sApvA_e;
        } else {
            ret->M = IvAp_e * sApvA_e;
        }
    } else {
        if (scale) {
            auto sIvAp_e(weight_matrix(IvAp_e, 0, true));
            ret->M = sIvAp_e * IvAp_e;
        } else {
            ret->M = std::move(IvAp_e);
        }
    }
    ret->dim[0] = std::move(dimI);
    ret->dim[1] = std::move(dimA);

    // Get weight vector from IvAp_e
    ret->weight.reference(sum_rows(IvAp_e));

    return ret;
}

static std::unique_ptr<WeightedSparse> compute_EvA(IceRegridder *sheet, bool scale, bool correctA, UrAE const &E, UrAE const &A)
{

    SparseSetT dimA, dimG, dimE;

    // ----- Get the Ur matrices (which determines our dense dimensions)

    // ApvG
    SparseTriplets<SparseMatrix> GvAp({&dimG, &dimA});
    GvAp.set_shape({sheet->nG(), A.nfull});
    A.GvAp(GvAp);

    // GvEp
    SparseTriplets<SparseMatrix> GvEp({&dimG, &dimE});
    GvEp.set_shape({sheet->nG(), E.nfull});
    E.GvAp(GvEp);


    // ----- Convert to Eigen and multiply
    auto GvAp_e(GvAp.to_eigen('.'));
    auto sGvAp_e(GvAp.eigen_scale_matrix(0));
    auto EpvG_e(GvEp.to_eigen('T'));

    // Unweighted matrix
    EigenSparseMatrix EpvAp_e(EpvG_e * sGvAp_e * GvAp_e);

    // ----- Apply final scaling, and convert back to sparse dimension
    std::unique_ptr<WeightedSparse> ret(new WeightedSparse);
    if (correctA) {

        // sApvA
        SparseTriplets<SparseMatrix> sApvA_t({&dimA, &dimA});
        sApvA_t.set_shape({A.nfull, A.nfull});
        A.sApvA(sApvA_t, in_sparse_fn(dimA));  // Filter to avoid adding more items to dimA
        auto sApvA_e(sApvA_t.to_eigen('.', false));       // Don't invert it...

        // ---- Obtain sEpvE
        SparseTriplets<SparseMatrix> sEpvE({&dimE, &dimE});
        sEpvE.set_shape({E.nfull, E.nfull});
        E.sApvA(sEpvE_t, std::bind(&SparseSetT::in_sparse, dimE, _1));    // Avoid adding more items to dimE

        if (scale) {
            auto sEvEp_e(sEpvE_t.to_eigen('.', true));        // Invert it...
            auto sEpvAp_e(weight_matrix(EpvAp_e, 0, true));

            ret->M = sEvEp_e * sEpvAp_e * EpvAp_e * sApvA_e;    // EvA
        } else {
            ret->M = EpvAp_e * sApvA_e;
        }

        // ----- Compute the final weight (diagonal) matrix
        auto wEvEp_e(sEpvE_t.to_eigen('.', false));
        auto wEpvAp_b(sum(EpvAp_e,0));    // _b = Blitz++
        EigenSparseMatrix weight_e(wEvEp_e * weight_matrix(wEpvAp_b));
        ret->weight.reference(wEpvAp_b);

    } else {    // ~correctA
        auto wEpvAp_b(sum(EpvAp_e,0));
        if (scale) {
            auto sEpvAp_e(scale_matrix(wEpvAp_b));

            ret->M = sEpvAp_e * EpvAp_e;
        } else {
            ret->M = std::move(EpvAp_e);
        }
        ret->dim[0] = std::move(dimE);
        ret->dim[1] = std::move(dimA);

        // ----- Compute the final weight (diagonal) matrix
        ret->weight.reference(wEpvAp_b);
    }

    return ret;

}

// ----------------------------------------------------------------
RegridMatrices::RegridMatrices(IceRegridder *sheet)
{
    printf("===== RegridMatrices Grid geometries:\n");
    printf("    nA = %d\n", sheet->gcm->nA());
    printf("    nhp = %d\n", sheet->gcm->nhp());
    printf("    nE = %d\n", sheet->gcm->nE());
    printf("    nI = %d\n", sheet->nI());
    printf("    nG = %d\n", sheet->nG());

    UrAE urA(sheet->gcm->nA(),
        std::bind(&IceRegridder::GvAp, sheet, _1),
        std::bind(&IceRegridder::sApvA, sheet, _1, _2));

    UrAE urE(sheet->gcm->nE(),
        std::bind(&IceRegridder::GvEp, sheet, _1),
        std::bind(&IceRegridder::sEpvE, sheet, _1, _2));

    // ------- AvI, IvA
    regrids.insert(make_pair("AvI", std::bind(&compute_AEvI, sheet, _1, _2, urA) ));
    regrids.insert(make_pair("IvA", std::bind(&compute_IvAE, sheet, _1, _2, urA) ));

    // ------- EvI, IvE
    regrids.insert(make_pair("EvI", std::bind(&compute_AEvI, sheet, _1, _2, urE) ));
    regrids.insert(make_pair("IvE", std::bind(&compute_IvAE, sheet, _1, _2, urE) ));

    // ------- EvA, AvE regrids.insert(make_pair("EvA", std::bind(&compute_EvA, this, _1, _2, urE, urA) ));
    regrids.insert(make_pair("EvA", std::bind(&compute_EvA, sheet, _1, _2, urE, urA) ));
    regrids.insert(make_pair("AvE", std::bind(&compute_EvA, sheet, _1, _2, urA, urE) ));

    // ----- Show what we have!
    printf("Available Regrids:");
    for (auto ii = regrids.begin(); ii != regrids.end(); ++ii) {
        printf(" %s", ii->first.c_str());
    }
    printf("\n");
}

}   // namespace icebin
