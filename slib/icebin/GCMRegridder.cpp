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
GCMRegridder_Standard::GCMRegridder_Standard() :
    mem_ice_regridders(mem_sheets_index)
{
    _ice_regridders = &mem_ice_regridders;
    mem_agridA.reset(new AbbrGrid);
    agridA = &*mem_agridA;
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
    mem_agridA.reset(new AbbrGrid(std::move(_agridA)));
    agridA = &*mem_agridA;
//    fgridA = std::move(_gridA);
//    agridA = AbbrGrid(*fgridA);
//  domainA = std::move(_domainA);
    this->_hcdefs = std::move(hcdefs);
    indexingHC = std::move(_indexingHC);
    correctA = _correctA;

    if (indexingHC.rank() != 2) (*icebin_error)(-1,
        "indexingHC has rank %d, it must have rank=2", indexingHC.rank());

    indexingE = derive_indexingE(agridA->indexing, indexingHC);
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
    agridA->clear();
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
    mem_agridA.reset(new AbbrGrid);
    agridA = &*mem_agridA;

    agridA->ncio(ncio, vname + ".agridA");
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


    indexingE = derive_indexingE(agridA->indexing, indexingHC);

}
// -------------------------------------------------------------


void GCMRegridder_Standard::filter_cellsA(std::function<bool(long)> const &keepA)
{

    // Now remove cells from the exgrids and gridIs that
    // do not interact with the cells we've kept in grid1.
    for (auto ice_regridder=ice_regridders().begin(); ice_regridder != ice_regridders().end(); ++ice_regridder) {
        (*ice_regridder)->filter_cellsA(keepA);
    }

    agridA->filter_cells(keepA);
}

void GCMRegridder_Standard::filter_cellsA(ibmisc::Domain const &domainA)
{
    filter_cellsA(std::bind(&ibmisc::in_domain,
        &domainA, &agridA->indexing, _1));
}

// ---------------------------------------------------------------------
// RegridMatrices const GCMRegridder_Standard::regrid_matrices(std::string const &sheet_name) const
//        ---> see RegridMatrices_Dynamic.cpp
// ---------------------------------------------------------------------


// ============================================================
// Special Debugging Functions

void check_nans(EigenSparseMatrixT const &M, std::string const &label)
{
    int nnan = 0;
    int ninf = 0;
    int nzero = 0;
    for (auto ii=begin(M); ii != end(M); ++ii) {
        if (std::isnan(ii->value())) ++nnan;
        if (std::isinf(ii->value())) ++ninf;
        if (ii->value() == 0) ++nzero;
    }
    printf("check_nans: Matrix %s has %d NaNs %d Infs %d zeros\n", label.c_str(), nnan, ninf, nzero);
}
void check_nans(blitz::Array<double,1> &M, std::string const &label)
{
    int nnan = 0;
    int ninf = 0;
    for (int i=0; i<M.extent(0); ++i) {
        if (std::isnan(M(i))) ++ nnan;
        if (std::isinf(M(i))) ++ninf;
    }
    printf("check_nans: Matrix %s has %d NaNs %d Infs\n", label.c_str(), nnan, ninf);
}


}   // namespace icebin
