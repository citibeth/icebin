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
#include <algorithm>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/cython.hpp>
#include <spsparse/SparseSet.hpp>
#include <icebin/IceRegridder.hpp>
#include <icebin/GCMCoupler.hpp>
#include <icebin/Grid.hpp>
#ifdef BUILD_MODELE
#include <icebin/modele/GCMCoupler_ModelE.hpp>
#endif
#include "icebin_cython.hpp"

using namespace ibmisc;
using namespace ibmisc::cython;
using namespace spsparse;
using namespace icebin;
using namespace icebin::modele;

namespace icebin {
namespace cython {

static double const nan = std::numeric_limits<double>::quiet_NaN();

void read_fgrid(
    std::unique_ptr<Grid> &fgridA,    // Return value goes here
    std::string const &fgridA_fname,
    std::string const &fgridA_vname)
{
    // Read fgridA
    NcIO ncio(fgridA_fname, netCDF::NcFile::read);
    fgridA.reset(new Grid);
    fgridA->ncio(ncio, fgridA_vname);
    ncio.close();
}

std::shared_ptr<GCMRegridder_Standard> new_GCMRegridder_Standard(
    Grid const &fgridA,
    std::vector<double> &hcdefs,
    bool _correctA)
{
    std::shared_ptr<GCMRegridder_Standard> cself(new GCMRegridder_Standard());

    // Put it together
    long nhc = hcdefs.size();
    cself->init(
        AbbrGrid(fgridA),
//        std::move(fgridA),
        std::move(hcdefs),
        Indexing({"A", "HC"}, {0,0}, {fgridA.ndata(), nhc}, {1,0}),
        _correctA);

    return cself;
}

std::shared_ptr<GCMRegridder> new_GCMRegridder_ModelE(
    std::string const &global_ecO,
    std::shared_ptr<GCMRegridder> const &gcmO)
{
#ifdef BUILD_MODELE
    std::shared_ptr<GCMRegridder_ModelE> gcmA(new GCMRegridder_ModelE(global_ecO, gcmO));
    return gcmA;
#else
    return std::shared_ptr<GCMRegridder_ModelE>();
#endif
}

void GCMRegridder_ModelE_set_focean(
    GCMRegridder *_gcmA,
    PyObject *foceanAOp_py,
    PyObject *foceanAOm_py)
{
#ifdef BUILD_MODELE
    auto gcmA(dynamic_cast<GCMRegridder_ModelE *>(_gcmA));

    // Check types and convert Numpy Arrays
    size_t nO = gcmA->gcmO->nA();
    auto _foceanAOp(np_to_blitz<double,1>(foceanAOp_py, "foceanAOp", {nO}));
    auto _foceanAOm(np_to_blitz<double,1>(foceanAOm_py, "foceanAOm", {nO}));

    // Copy values from Python memory to C++ memory
    gcmA->foceanAOp = _foceanAOp;
    gcmA->foceanAOm = _foceanAOm;
#endif
}

PyObject *GCMRegridder_wA(
    GCMRegridder *gcm_regridder,
    std::string const &sheet_name,
    bool native,    // native vs. projected grid
    double fill)
{
    PyObject *wA_py = ibmisc::cython::new_pyarray<double,1>(
        std::array<int,1>{gcm_regridder->agridA.dim.sparse_extent()});
    auto wA(np_to_blitz<double,1>(wA_py, "wA", {-1}));
    wA = fill;

    gcm_regridder->wA(
        accum::blitz_existing(wA, DuplicatePolicy::REPLACE), sheet_name, native);

    return wA_py;
}


void GCMRegridder_add_sheet(GCMRegridder *cself,
    Grid const &fgridA,
    std::string const &name,
    std::string const &gridI_fname, std::string const &gridI_vname,
    std::string const &exgrid_fname, std::string const &exgrid_vname,
    std::string const &sinterp_style)
{
    NcIO ncio_I(gridI_fname, netCDF::NcFile::read);
    std::unique_ptr<Grid> fgridI(new Grid);
    fgridI->ncio(ncio_I, "grid");
    ncio_I.close();

    NcIO ncio_exgrid(exgrid_fname, netCDF::NcFile::read);
    std::unique_ptr<Grid> fexgrid(new Grid);
    fexgrid->ncio(ncio_exgrid, exgrid_vname);
    ncio_exgrid.close();

    auto interp_style(parse_enum<InterpStyle>(sinterp_style));

    auto sheet(new_ice_regridder(fgridI->parameterization));
    sheet->init(
        name, cself->agridA, &fgridA,
        AbbrGrid(*fgridI), ExchangeGrid(*fexgrid),
        interp_style);

    dynamic_cast<GCMRegridder_Standard *>(cself)
        ->add_sheet(std::move(sheet));
}

/** Computes yy = M xx.  yy is allocated, not necessarily set. */
void coo_matvec(PyObject *yy_py, PyObject *xx_py, bool ignore_nan,
    size_t M_nrow, size_t M_ncol, PyObject *M_row_py, PyObject *M_col_py, PyObject *M_data_py)
{
    auto xx(np_to_blitz<double,1>(xx_py, "xx", {M_ncol}));
    auto yy(np_to_blitz<double,1>(yy_py, "yy", {M_nrow}));
    auto M_row(np_to_blitz<int,1>(M_row_py, "M_row_py", {-1}));
    auto M_col(np_to_blitz<int,1>(M_col_py, "M_col_py", {-1}));
    auto M_data(np_to_blitz<double,1>(M_data_py, "M_data_py", {-1}));

    // Keep track of which items we've written to.
    std::vector<bool> written(M_nrow, false);

    // Do the multiplication, and we're done!
    int nnz = M_data.size();
    for (int n=0; n<nnz; ++n) {
        size_t row = M_row(n);
        size_t col = M_col(n);
        double data = M_data(n);

        // Ignore NaN in input vector
        if (ignore_nan && std::isnan(xx(col))) continue;

        // Just do Snowdrift-style "REPLACE".  "MERGE" was never used.
        double old_yy;
        if (written[row]) {
            old_yy = yy(row);
        } else {
            old_yy = 0;
            written[row] = true;
        }
        yy(row) = old_yy + data * xx(col);
    }

}


extern linear::Weighted *RegridMatrices_matrix(RegridMatrices *cself, std::string const &spec_name)
{
    return cself->matrix(spec_name).release();
}


// ------------------------------------------------------------
PyObject *Hntr_regrid(Hntr const *hntr, PyObject *WTA_py, PyObject *A_py, bool mean_polar)
{
    // Get Fortran-style (base index = 1) arrays out of this.
    auto WTA(np_to_blitz<double,1>(WTA_py, "WTA", {hntr->Agrid.spec.size()}, blitz::fortranArray));
    auto A(np_to_blitz<double,1>(A_py, "WTA", {hntr->Agrid.spec.size()}, blitz::fortranArray));
    PyObject *B_py(ibmisc::cython::new_pyarray<double,1>(make_array(hntr->Bgrid.spec.size())));
    auto B(np_to_blitz<double,1>(B_py, "WTA", {-1}, blitz::fortranArray));

    hntr->regrid(WTA, A, B, mean_polar);

    return B_py;
}

RegridMatrices *new_regrid_matrices(
    GCMRegridder const *gcm,
    std::string const &sheet_name,
    PyObject *elevmaskI_py,
    // --------- Params
    bool scale,
    bool correctA,
    double sigma_x,
    double sigma_y,
    double sigma_z,
    bool conserve)
{


    auto sheet_index = gcm->ice_regridders().index.at(sheet_name);
    IceRegridder *ice_regridder = &*gcm->ice_regridders()[sheet_index];
    auto elevmaskI(np_to_blitz<double,1>(elevmaskI_py, "elevmaskI", {ice_regridder->nI()}));

    return gcm->regrid_matrices(
        sheet_index, elevmaskI,
        RegridParams(scale, correctA, {sigma_x, sigma_y, sigma_z})).release();
}

std::string to_string(PyObject *str, std::string const &vname)
{
    if (!PyUnicode_Check(str)) (*icebin_error)(-1,
        "Variable %s must be str", vname.c_str());

    Py_ssize_t size;
    char *buf = PyUnicode_AsUTF8AndSize(str, &size);
    return std::string(buf, size);
}

#if 0
void update_topo(
    // ====== INPUT parameters
    GCMRegridder *_gcmA,
    std::string const &topoO_fname,    // Name of Ocean-based TOPO file (aka Gary)
    PyObject *elevmask_sigmas_py,    // {'greenland' : (elevmaskI<1>, maskI<1>, (sigma_x,signa_y,sigma_z)), ...}
    bool initial_timestep,    // true if this is the first (initialization) timestep
    std::string const &segments,    // string, eg: "legacy,sealand,ec"
    std::string const &primary_segment,    // [('name', base), ...]
    // ===== OUTPUT parameters (variables come from GCMCoupler); must be pre-allocated
    PyObject *fhc_py,         // blitz::Array<double,3> fhc;
    PyObject *underice_py,    // blitz::Array<int,3> underice;
    PyObject *elevE_py,       // blitz::Array<double,3> elevE;
    // i,j arrays on Atmosphere grid
    PyObject *focean_py,      // blitz::Array<double,2> focean;
    PyObject *flake_py,       // blitz::Array<double,2> flake;
    PyObject *fgrnd_py,       // blitz::Array<double,2> fgrnd;    // Alt: fearth0
    PyObject *fgice_py,       // blitz::Array<double,2> fgice;    // Alt: flice
    PyObject *zatmo_py,       // blitz::Array<double,2> zatmo;      // i,j

    PyObject *foceanOm0_py)   // blitz::Array<double,1> foceanOm0,
{
#ifdef BUILD_MODELE
    auto gcmA(dynamic_cast<GCMRegridder_ModelE *>(_gcmA));

    // --------------------------------------------------------
    // Convert elevmask and sigmas to C++ Data Structures
    std::map<std::string, std::pair<ElevMask<1>, std::array<double,3>>> tdata;
    if (!PyDict_Check(elevmask_sigmas_py)) (*icebin_error)(-1,
        "elevmask_sigmas must be a Python dict.");
    PyObject *sheet_py, *value;
    Py_ssize_t ppos = 0;
    while (PyDict_Next(elevmask_sigmas_py, &ppos, &sheet_py, &value)) {
        // Get nI
        std::string sheet(to_string(sheet_py, "sheet_py"));
        IceRegridder *icer = &*gcmA->ice_regridders().at(sheet);
        int nI = icer->nI();

        // Parse the main tuple
        PyObject *elevmaskI_py, *maskI_py, *sigma_py;
        PyArg_ParseTuple(value, "OOO", &elevmaskI_py, &maskI_py, &sigma_py);
        auto elevmaskI(np_to_blitz<double,1>(elevmaskI_py, "elevmaskI", {nI}));
        auto maskI(np_to_blitz<char,1>(maskI_py, "maskI", {nI}));

        // Further parse sigma_py
        std::array<double,3> sigma;
        PyArg_ParseTuple(sigma_py, "ddd", &sigma[0], &sigma[1], &sigma[2]);

        // Add to temporary C++ Data Structures
        tdata.insert(std::make_pair(sheet,
            std::make_pair(ElevMask<1>(elevmaskI, maskI), sigma)));
    }

    // Convert temporary to permanent C++ data structure
    std::vector<ElevMask<1>> elevmasks;
    std::vector<std::array<double,3>> sigmas;
    auto &index(gcmA->ice_regridders().index);
    for (size_t i=0; i<index.size(); ++i) {
        auto const &rec(tdata.find(index[i]));
        elevmasks.push_back(rec->second.first);
        sigmas.push_back(rec->second.second);
    }

    // --------------------------------------------------------
    // Convert segments to C++ Data Structure
    std::vector<HCSegmentData> hc_segments(parse_hc_segments(segments));
    auto const &ec(get_segment(hc_segments, "ec"));
    auto nhc_ice(gcmA->nhc());
    int nhc_gcm = ec.base + nhc_ice;

    // --------------------------------------------------------
    // Convert arrays to C++ Data Structures
    GridSpec_LonLat *specA = dynamic_cast<GridSpec_LonLat *>(&*gcmA->agridA.spec);
    int nj = specA->nlat();
    int ni = specA->nlon();

    icebin::modele::Topos toposA;

    // Convert 3D on the Elevation Grid
    toposA.fhc.reference(np_to_blitz<double,3>(fhc_py, "fhc", {nhc_gcm, nj, ni}));
    toposA.underice.reference(np_to_blitz<int,3>(underice_py, "underice", {nhc_gcm, nj, ni}));
    toposA.elevE.reference(np_to_blitz<double,3>(elevE_py, "elevE", {nhc_gcm, nj, ni}));

    // Convert 2D on the Atmosphere Grid
    toposA.focean.reference(np_to_blitz<double,2>(focean_py, "focean", {nj, ni}));
    toposA.flake.reference(np_to_blitz<double,2>(flake_py, "flake", {nj, ni}));
    toposA.fgrnd.reference(np_to_blitz<double,2>(fgrnd_py, "fgrnd", {nj, ni}));
    toposA.fgice.reference(np_to_blitz<double,2>(fgice_py, "fgice", {nj, ni}));
    toposA.zatmo.reference(np_to_blitz<double,2>(zatmo_py, "zatmo", {nj, ni}));

    // Convert 2D on the Ocean Grid
    GCMRegridder *gcmO = &*gcmA->gcmO;
    GridSpec_LonLat *specO = dynamic_cast<GridSpec_LonLat *>(&*gcmO->agridA.spec);
    int njO = specO->nlat();
    int niO = specO->nlon();
    auto foceanOm0(np_to_blitz<double,2>(foceanOm0_py, "foceanOm0", {njO, niO}));

    // --------------------------------------------------------
    // Call C++ Function
    icebin::modele::update_topo(gcmA, topoO_fname, elevmasks, sigmas,
        initial_timestep, hc_segments, primary_segment, toposA, foceanOm0);
#endif
}
#endif    // if 0

}}
