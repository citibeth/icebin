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

    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        ice_couplers[sheetix]->set_start_time(time_base, time_start_s);
    }

}
// ------------------------------------------------------------

GCMCouplerOutput GCMCoupler::couple(
// Simulation time [s]
double time_s,
ArraySparseParallelVectorsE const &gcm_ovalsE,
bool do_run)
{
    // Initialize output
    GCMCouplerOutput out;
    for (int i=0; i<GCMI::COUNT; ++i) {
        out.gcm_ivals[i].nvar = gcm_inputs[i].size();
    }

    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        auto &ice_coupler(ice_couplers[sheetix]);
        ice_coupler.couple(time_s, gcm_ovalsE, out, do_run);
    }
    return out;
}


}       // namespace
