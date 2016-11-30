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

#include <mpi.h>        // Intel MPI wants to be first
#include <ibmisc/netcdf.hpp>
#include <ibmisc/memory.hpp>
#include <ibmisc/blitz.hpp>
#include <icebin/GCMCoupler.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/contracts/contracts.hpp>
#include <spsparse/multiply_sparse.hpp>
#include <spsparse/sort.hpp>

#ifdef USE_PISM
#include <icebin/pism/IceModel_PISM.hpp>
#endif

using namespace spsparse;
using namespace ibmisc;
using namespace netCDF;

namespace icebin {

// ==========================================================

#if 0
/** Converts a VectorSparseParallelVectors structure to an
    ArraySparseParallelVectors (used as arg to functions). */
ArraySparseParallelVectors vector_to_array(VectorSparseParallelVectors vecs)
{
    ArraySparseParallelVectors ret;
    ret.index.reference(to_blitz(vecs.index));

    blitz::TinyVector<int,1> shape(0);
    blitz::TinyVector<int,1> strides(0);

    shape[0] = vec.size();
    strides[0] = vecs.nvar;     /* Blitz++ strides in sizeof(T) units */

    for (int i=0; i<vecs.nvar; ++i) {
        ret.vals.push_back(&vec[i], shape, strides,
            blitz::neverDeleteData)
    }
    return ret;
}
#endif
// ==========================================================
/** @param nc The IceBin configuration file */
void GCMCoupler::ncread(
    std::string const &_fname,
    std::string const &_vname,
    ibmisc::Domain<int> &&_domainA)
{
    printf("BEGIN GCMCoupler::read_from_netcdf() %s\n", vname.c_str());

    fname = _fname;
    vname = _vname;
    domainA = std::move(_domainA);

    NcIO ncio(fname, NcFile::read);

    // Load the MatrixMaker (filtering by our domain, of course)
    // Also load the ice sheets
    regridder.ncio(ncio, vname);
    regridder.filter_cellsA(domainA);

    // The root node needs to make the full regridding matrices, not just those
    // for its own domain.  Therefore, create a second MatrixMaker for it.
    if (am_i_root()) {
        regridder_full.reset(new GCMRegridder);
        regridder_full->ncio(ncio, vname);
    }

    // Read gcm_out_file, an optional variable telling the GCM-specific
    // part of IceBin to write out exactly what it sees coming from the GCM
    // (so it can be replayed later with desm)
    auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});

    gcm_out_file = "";
    get_or_put_att(info_v, ncio.rw, "gcm_out_file", gcm_out_file, false);
    if (gcm_out_file.length() > 0) {
        gcm_out_file = boost::filesystem::absolute(
            boost::filesystem::path(gcm_out_file),
            gcm_params.run_dir).string();
    }

    gcm_in_file = "";
    get_or_put_att(info_v, ncio.rw, "gcm_in_file", gcm_in_file, false);
    if (gcm_in_file.length() > 0) {
        gcm_in_file = boost::filesystem::absolute(
            boost::filesystem::path(gcm_in_file),
            gcm_params.run_dir).string();
    }



#if 1
    std::cout << "========= GCM Constants" << std::endl;
    std::cout << gcm_constants;
    std::cout << "========= GCM Outputs" << std::endl;
    std::cout << gcm_outputs;
    std::cout << "========= GCM Inputs" << std::endl;
    std::cout << gcm_inputs;


#endif

    ice_couplers.clear();
    for (size_t i=0; i < regridder.sheets.size(); ++i) {
        IceRegridder *sheet = &*regridder.sheets[i];
        std::string vname_sheet = vname + "." + sheet->name();

        // Create an IceCoupler corresponding to this IceSheet.
        std::unique_ptr<IceCoupler> ice_coupler(new_ice_coupler(ncio, vname_sheet, this, sheet));

        contracts::setup(*this, *ice_coupler);    // where does this go w.r.t ncread() and upate?
        ice_coupler->update_ice_sheet(ncio, vname_sheet);

        ice_couplers.push_back(std::move(ice_coupler));
    }

    ncio.close();
    printf("END GCMCoupler::read_from_netcdf()\n");
}



void GCMCoupler::set_start_time(
    ibmisc::time::tm const &time_base,
    double time_start_s)
{
    gcm_params.set_start_time(time_base, time_start_s);
    last_time_s = time_start_s;

    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        ice_couplers[sheetix]->set_start_time(time_base, time_start_s);
    }

}
// ------------------------------------------------------------

GCMCouplerOutput GCMCoupler::couple(
// Simulation time [s]
double time_s,
std::array<int,3> const &yymmdd, // Date that time_s lies on
ArraySparseParallelVectorsE const &gcm_ovalsE,
bool do_run)
{
    std::string sdate = (boost::format
        ("%04d%02d%02d") % yymmdd[0] % yymmdd[1] % yymmdd[2]).str();
    std::string log_dir = "icebin";

    if (gcm_params.icebin_logging) {
        std::string fname = "gcm-out-" + sdate;
        NcIO ncio(fname, 'w');
        ncio_gcm_output(ncio, gcm_ovalsE, {last_time_s, time_s},
            gcm_params.time_units, "");
        ncio();
    }

    // Initialize output
    GCMCouplerOutput out(gcm_inputs[i].size());

    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        auto &ice_coupler(ice_couplers[sheetix]);
        ice_coupler.couple(time_s, gcm_ovalsE, out, do_run);
    }


    if (gcm_params.icebin_logging) {
        std::string fname = "gcm-out-" + sdate;
        NcIO ncio(fname, 'r');
        ncio_gcm_input(ncio, gcm_ovalsE, {last_time_s, time_s},
            gcm_params.time_units, "");
        ncio();
    }


    return out;
}
// ------------------------------------------------------------
static ncwrite_dense(
    NcFile *nc,
    ArraySparseParallelVectors *vecs,
    VarSet *contract,
    ibmisc::IndexingBase const *indexing,
    std::string const &vname_base)
{
    // im,jm,ihc  0-based
    blitz::Array<double,1> denseE(indexing->extent());    // All dims

    // NetCDF needs dimensions in stride-descending order
    std::vector<size_t> startp
    std::vector<size_t> extents;
    for (int dimi : indexing->indices) {
        startp.push_back(0);
        extents.push_back(indexing->extents[dimi]);
    }

    // Go through each variable...
    for (unsigned int ivar=0; ivar < vals.values.size(); ++ivar) {
        // Fill our dense var
        denseE = nan;
        for (unsigned int i=0; i<vals.index.size(); ++i)
            dense(vals.index[i]) = vals.values[ivar][i];

        // Store in the netCDF variable
        NcVar ncvar(nc.getVar(vname_base + (*contract)[ivar].name));
        ncvar.putVar(startp, countp, denseE.data());
    }

}

/** Densifies the sparse vectors, and then writes them out to netCDF */
void ncio_dense(
    NcIO &ncio,
    ArraySparseParallelVectors &vec,    
    VarSet &contract,
    ibmisc::IndexingBase const &indexing,
    std::string const &vname_base = "")
{
    if (ncio.rw != 'w') (*icebin_error)(-1,
        "ncio_dense(ArraySparseParallelVectors) only writes, no read.");

    // Create the variables
    NcDimSpec dim_spec({"time"}, {-1});
    append(dim_spec, indexingAE);
    contract.ncdefine(ncio, dim_spec.to_dims(), vname_base);

    ncio += std::bind(&ncwrite_dense, ncio.nc,
        &vec, &contract, &indexing,
        vname_base);
}
// ------------------------------------------------------------
template<ValT>
inline append(std::vector<ValT> &out, std::vector<ValT> const &inn)
{
    for (auto ii=inn.begin(); ii != inn.end(); ++ii)
        out.push_back(*ii);
}
GCMCouplerOutput concatenate(std::vector<GCMCouplerOutput> const &outs)
{
    GCMCouplerOutput all_out;

    if (outs.size() == 0) return all_out;

    for (GridAE iAE=0; iAE<GridAE::count; ++iAE) {
        append(all_out.gcm_ivals[iAE].index, out.gcm_ivals[iAE].index);
        append(all_out.gcm_ivals[iAE].vals, out.gcm_ivals[iAE].vals);
    }

    all_out = outs[0];
    for (size_t i=1; i<outs.size(); ++i) {
        GCMCouplerOutput &out(outs[i]);

        copy(all_out.E1vE0, E1vE0);
        copy(all_out.AvE1, AvE1);
        copy(all_out.wAvE1, wAvE1);
        copy(all_out.elevE1, elevE1);
    }

    return all_out;
}
// ------------------------------------------------------------
/** Top-level ncio() to log output from coupler. (coupler->GCM) */
void GCMCoupler::ncio_gcm_input(NcIO &ncio,
    GCMCouplerOutput &out,        // All MPI ranks included here
    std::array<double,2> const &timespan,    // timespan[1] = current time_s
    std::string time_units,
    std::string const &vname_base)
{
    ncio_timespan(ncio, timespan, time_units, vname_base);

    for (GridAE iAE=0; iAE<GridAE::count; ++iAE) {
        ncio_dense(ncio, out.gcm_ivals[iAE], gcm_inputs[iAE],
            regridder.indexing(iAE), vname_base);
    }

    ncio_spsparse(ncio, out.E1vE0, false, vname_base+"Ev1Ev0");
    ncio_spsparse(ncio, out.AvE1, false, vname_base+"AvE1");
    ncio_spsparse(ncio, out.wAvE1, false, vname_base+"wAvE1");
    ncio_spsparse(ncio, out.elevE1, false, vname_base+"elevE1");
}

/** Top-level ncio() to log input to coupler. (GCM->coupler) */
void GCMCoupler::ncio_gcm_output(NcIO &ncio,
    ArraySparseParallelVectors &gcm_ovalsE,
    std::array<double,2> const &timespan,    // timespan[1] = current time_s
    std::string time_units,
    std::string const &vname_base)
{
    ncio_timespan(ncio, timespan, time_units, vname_base);

    ncio_dense(ncio, gcm_ovalsE, gcm_outputsE,
        regridder.indexingE, vname_base);
}
// ------------------------------------------------------------

}       // namespace
