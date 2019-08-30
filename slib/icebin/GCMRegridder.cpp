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
#include <spsparse/netcdf.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/Grid.hpp>

using namespace std;
using namespace netCDF;
using namespace ibmisc;
using namespace std::placeholders;  // for _1, _2, _3...
using namespace spsparse;

namespace icebin {

GCMRegridder::~GCMRegridder() {}


// ========================================================
/** Determines correct indexing for the E grid.
Returns: im,jm,...,iHC*/
ibmisc::Indexing derive_indexingE(
ibmisc::Indexing const &indexingA,      // im,jm,...
ibmisc::Indexing const &indexingHC)    // iA,iHC
{
    std::vector<IndexingData> dims(indexingA.dims());
    std::vector<int> indices(indexingA.indices());   // Index IDs sorted by descending stride. {0,1,...} for row-major, reversed for col-major

    dims.push_back(indexingHC[1]);

    if (indexingHC.indices()[0] == 0) {    // HC is row major so we are too
        indices.push_back(indexingA.rank());
    } else {
        for (auto &ii : indices) ++ii;
        indices.push_back(0);
    }

    return ibmisc::Indexing(std::move(dims), std::move(indices));
}
// -------------------------------------------------------------
void GCMRegridder_Standard::init(
    AbbrGrid &&_agridA,
//    std::unique_ptr<Grid> &&_gridA,
//  ibmisc::Domain<int> &&_domainA,     // Tells us which cells in gridA to keep...
    std::vector<double> &&hcdefs,  // [nhc]
    ibmisc::Indexing &&_indexingHC,
    bool _correctA)
{
    agridA = std::move(_agridA);
//    fgridA = std::move(_gridA);
//    agridA = AbbrGrid(*fgridA);
//  domainA = std::move(_domainA);
    this->_hcdefs = std::move(hcdefs);
    indexingHC = std::move(_indexingHC);
    correctA = _correctA;

    if (indexingHC.rank() != 2) (*icebin_error)(-1,
        "indexingHC has rank %d, it must have rank=2", indexingHC.rank());

    indexingE = derive_indexingE(agridA.indexing, indexingHC);
}
// -------------------------------------------------------------

// -------------------------------------------------------------
// ==============================================================
void GCMRegridder::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    (*icebin_error)(-1, "Not implemented");
}
// ==============================================================
void GCMRegridder_Standard::clear()
{
    agridA.clear();
    _hcdefs.clear();
    ice_regridders().index.clear();
    ice_regridders().clear();
}
// -------------------------------------------------------------
void GCMRegridder_Standard::ncio(NcIO &ncio, std::string const &vname)
{
    auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});

    if (ncio.rw == 'r') {
        clear();
    }

    // Read/Write gridA and other global stuff
    agridA.ncio(ncio, vname + ".agridA");
    indexingHC.ncio(ncio, vname + ".indexingHC");
    ncio_vector(ncio, _hcdefs, true, vname + ".hcdefs", "double",
        get_or_add_dims(ncio, {vname + ".nhc"}, {_hcdefs.size()} ));
//  domainA.ncio(ncio, ncInt, vname + ".domainA");
    get_or_put_att(info_v, ncio.rw, "correctA", &correctA, 1);

    // Read/write list of sheet names
    std::vector<std::string> sheet_names;
    if (ncio.rw == 'w') {
        // Copy sheet names to a std::vector
        for (auto ii=ice_regridders().index.begin(); ii != ice_regridders().index.end(); ++ii)
            sheet_names.push_back(*ii);
    }
    get_or_put_att(info_v, ncio.rw, "sheets", "string", sheet_names);

    // Read/write the regridder for each ice sheet
    if (ncio.rw == 'r') {   // Instantiate
        for (auto sheet_name = sheet_names.begin(); sheet_name != sheet_names.end(); ++sheet_name) {
            std::string vn(vname + "." + *sheet_name);
            add_sheet(*sheet_name, new_ice_regridder(ncio, vn));
        }
    }
    for (auto ice_regridder=ice_regridders().begin(); ice_regridder != ice_regridders().end(); ++ice_regridder) {
        (*ice_regridder)->ncio(ncio, vname + "." + (*ice_regridder)->name());
    }


    indexingE = derive_indexingE(agridA.indexing, indexingHC);

}
// -------------------------------------------------------------


void GCMRegridder_Standard::filter_cellsA(std::function<bool(long)> const &keepA)
{

    // Now remove cells from the exgrids and gridIs that
    // do not interact with the cells we've kept in grid1.
    for (auto ice_regridder=ice_regridders().begin(); ice_regridder != ice_regridders().end(); ++ice_regridder) {
        (*ice_regridder)->filter_cellsA(keepA);
    }

    agridA.filter_cells(keepA);
}

void GCMRegridder_Standard::filter_cellsA(ibmisc::Domain const &domainA)
{
    filter_cellsA(std::bind(&ibmisc::in_domain,
        &domainA, &agridA.indexing, _1));
}

// ---------------------------------------------------------------------
// RegridMatrices const GCMRegridder_Standard::regrid_matrices(std::string const &sheet_name) const
//        ---> see RegridMatrices_Dynamic.cpp
// ---------------------------------------------------------------------
#if 0
/** Default implementation. */
linear::Weighted_Tuple GCMRegridder_Standard::global_unscaled_AvE(
    std::vector<blitz::Array<double,1>> const &emI_lands,
    std::vector<blitz::Array<double,1>> const &emI_ices,
    RegridParams const &params) const
{
    (*icebin_error)(-1, "Generic GCMRegridder::global_unscaled_AvE() has not been tested; but it should be close to correct, if needed.");
    linear::WeightedTuple AvE_g;    // _g = global


    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        auto &ice_sheet(*ice_regridders()[sheet_index]);
        auto rm(regrid_matrices(sheet_index, emI_ices[sheet_index], params));

        SparseSetT dimA1;
        SparseSetT dimE;
        // ---- Update AvE matrix and weights (global for all ice sheets)
        // Adds to dimA1 and dimE
        auto AvE(rm->matrix_d("AvE", {&dimA1, &dimE}, regrid_params));

        // Accumulate per-ice sheet matrix and weights into total.
        spcopy(
            accum::to_sparse({AvE->dims[0]}),
            accum::ref(AvE_g.wM),
            AvE.wM);

        spcopy(
            accum::to_sparse(AvE->dims,
            accum::ref(AvE_g.M)),
            *AvE.M);

        spcopy(
            accum::to_sparse({AvE->dims[1]}),
            accum::ref(AvE_g.Mw),
            AvE.Mw);
    }    // Flush accumulators at destruction time

    ret->set_shape(std::array<long,2>{nA(), nE()});
    return AvE;
}
#endif
// ------------------------------------------------------------
linear::Weighted_Tuple GCMRegridder_Standard::global_unscaled_E1vE0(
    std::vector<linear::Weighted_Eigen *> const &E1vIs_unscaled, // State var set in IceCoupler::couple(); _nc = no correctA (RegridParam) UNSCALED matrix
    std::vector<EigenSparseMatrixT *> const &IvE0s, // State var set in IceCoupler::couple()  SCALED matrix
    std::vector<SparseSetT *> const &dimE0s) const    // dimE0 accompanying IvE0
{
    linear::Weighted_Tuple E1vE0_g;    // _g = global

    // ---------- Compute each E1vE0 and merge...
    for (size_t sheet_index=0; sheet_index < ice_regridders().size(); ++sheet_index) {
        auto &E1vI_unscaled(*E1vIs_unscaled[sheet_index]);
        EigenSparseMatrixT &IvE0(*IvE0s[sheet_index]);
        auto &dimE0(*dimE0s[sheet_index]);

        // Don't do this on the first round, since we don't yet have an IvE0
        EigenSparseMatrixT E1vE0(*E1vI_unscaled.M * IvE0);    // UNSCALED
        spcopy(
            accum::to_sparse(make_array(E1vI_unscaled.dims[0]),
            accum::ref(E1vE0_g.wM)),    // Output
            E1vI_unscaled.wM);   // Input

        spcopy(
            accum::to_sparse(make_array(E1vI_unscaled.dims[0], &dimE0),
            accum::ref(E1vE0_g.M)),
            E1vE0);

#if 0    // Not needed
        spcopy(
            accum::to_sparse(make_array(&dimE0),
            accum::ref(E1vE0_g.Mw)),
            IvE0.Mw);
#endif
    }    // Flush accumulators

    E1vE0_g.set_shape(std::array<long,2>{nE(), nE()});
    return E1vE0_g;
}


}   // namespace icebin
