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
#include <functional>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <ibmisc/netcdf.hpp>
#include <ibmisc/ncfile.hpp>
#include <ibmisc/memory.hpp>
#include <ibmisc/blitz.hpp>
#include <ibmisc/datetime.hpp>
#include <icebin/GCMCoupler.hpp>
#include <icebin/GCMRegridder.hpp>
#include <icebin/contracts/contracts.hpp>
#include <spsparse/netcdf.hpp>

#ifdef USE_PISM
#include <icebin/pism/IceCoupler_PISM.hpp>
#endif

using namespace spsparse;
using namespace ibmisc;
using namespace netCDF;

namespace icebin {

static double const nan = std::numeric_limits<double>::quiet_NaN();

// ==========================================================
static std::vector<HCSegmentData> parse_hc_segments(std::string const &str)
{
    std::vector<HCSegmentData> ret;

    // TDOO: Parse for real later
    std::vector<std::string> segment_names;
    boost::algorithm::split(segment_names, str, boost::is_any_of(","));
    int base = 0;
    for (auto &seg : segment_names) {
        boost::algorithm::trim(seg);
        int len;
        if (seg == "legacy") len = 1;
        else if (seg == "sealand") len = 2;
        else len = -1;

        ret.push_back(HCSegmentData(seg, base, len));
        base += len;
    }
    return ret;
}


/** @param nc The IceBin configuration file */
void GCMCoupler::ncread(
    std::string const &grid_fname,        // comes from this->gcm_params
    std::string const &config_fname,        // comes from this->gcm_params
    std::string const &vname)        // comes from this->gcm_params
{
    printf("BEGIN GCMCoupler::ncread(%s)\n", grid_fname.c_str()); fflush(stdout);

    // Load the MatrixMaker (filtering by our domain, of course)
    // Also load the ice sheets
    {
        NcIO ncio(grid_fname, NcFile::read);
        gcm_regridder.ncio(ncio, vname);
    }

    std::cout << "========= GCM Constants" << std::endl;
    std::cout << gcm_constants;
    std::cout << "========= GCM Outputs" << std::endl;
    std::cout << gcm_outputsE;
    std::cout << "========= GCM InputsA" << std::endl;
    std::cout << gcm_inputsAE[GridAE::A];
    std::cout << "========= GCM InputsA" << std::endl;
    std::cout << gcm_inputsAE[GridAE::E];

    NcIO ncio(config_fname, NcFile::read);

    // Read m.segments
    std::string segments;
    auto info_v = get_or_add_var(ncio, vname + ".info", "int64", {});
    get_or_put_att(info_v, ncio.rw, "segments", segments);
    gcm_params.hc_segments = parse_hc_segments(segments);
    gcm_params.icebin_base_hc = gcm_params.ec_segment().base;

    ice_couplers.clear();
    for (size_t i=0; i < gcm_regridder.ice_regridders.size(); ++i) {
        IceRegridder *ice_regridder = &*gcm_regridder.ice_regridders[i];
        std::string vname_sheet = vname + "." + ice_regridder->name();

        // Create an IceCoupler corresponding to this IceSheet.
        std::unique_ptr<IceCoupler> ice_coupler(new_ice_coupler(ncio, vname_sheet, this, ice_regridder));

        ice_couplers.push_back(std::move(ice_coupler));
    }

    printf("END GCMCoupler::ncread(%s)\n", grid_fname.c_str()); fflush(stdout);
}



void GCMCoupler::cold_start(
    ibmisc::Datetime _time_base,
    double _time_start_s)
{
    time_base = _time_base;
    time_start_s = _time_start_s;
    time_unit = TimeUnit(&cal365, time_base, TimeUnit::SECOND);
    last_time_s = time_start_s;

    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        auto &ice_coupler(ice_couplers[sheetix]);

        // NOTE: Contract sizes are only needed at runtime, not allocate() time.
        contracts::setup(*this, *ice_coupler);    // where does this go w.r.t ncread() and upate?
        ice_coupler->cold_start(time_base, time_start_s);
    }

}
// ------------------------------------------------------------

GCMInput GCMCoupler::couple(
double time_s,        // Simulation time [s]
VectorMultivec const &gcm_ovalsE,
bool run_ice,
bool am_i_root)
{
    if (!am_i_root) {
        GCMInput out({0,0});
        for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
            auto &ice_coupler(ice_couplers[sheetix]);
            ice_coupler->couple(time_s, gcm_ovalsE, out, run_ice, am_i_root);
        }
        return out;
    }

    std::array<double,2> timespan{last_time_s, time_s};

    // Figure out our calendar day to format filenames
    ibmisc::Datetime dt(time_unit.to_datetime(time_s));

    std::string sdate = (boost::format
        ("%04d%02d%02d") % dt[0] % dt[1] % dt[2]).str();
    std::string log_dir = "icebin";

    if (gcm_params.icebin_logging) {
        std::string fname = "gcm-out-" + sdate;
        NcIO ncio(fname, 'w');
        ncio_gcm_output(ncio, gcm_ovalsE, timespan,
            time_unit.to_cf(), "");
        ncio();
    }

    // Initialize output
    GCMInput out({gcm_inputsAE[0].size(), gcm_inputsAE[1].size()});

    for (size_t sheetix=0; sheetix < ice_couplers.size(); ++sheetix) {
        auto &ice_coupler(ice_couplers[sheetix]);
        ice_coupler->couple(time_s, gcm_ovalsE, out, run_ice, am_i_root);
    }


    if (gcm_params.icebin_logging) {
        std::string fname = "gcm-in-" + sdate;
        NcIO ncio(fname, 'r');
        ncio_gcm_input(ncio, out, timespan,
            time_unit.to_cf(), "");
        ncio();
    }

    return out;
}
// ------------------------------------------------------------
static void ncwrite_dense(
    NcGroup * const nc,
    VectorMultivec const *vecs,
    VarSet const *contract,
    ibmisc::Indexing const *indexing,
    std::string const &vname_base)
{
    // im,jm,ihc  0-based
    blitz::Array<double,1> denseE(indexing->extent());    // All dims

    // NetCDF needs dimensions in stride-descending order
    std::vector<size_t> startp;
    std::vector<size_t> countp;
    for (int dimi : indexing->indices()) {
        startp.push_back(0);
        countp.push_back(indexing[dimi].extent());
    }

    // Go through each variable...
    unsigned int nvar = vecs->size();
    for (unsigned int ivar=0; ivar < nvar; ++ivar) {
        // Fill our dense var
        denseE = nan;
        for (unsigned int i=0; i<vecs->index.size(); ++i)
            denseE(vecs->index[i]) = vecs->vals[i*nvar + ivar];

        // Store in the netCDF variable
        NcVar ncvar(nc->getVar(vname_base + (*contract)[ivar].name));
        ncvar.putVar(startp, countp, denseE.data());
    }

}

/** Densifies the sparse vectors, and then writes them out to netCDF */
static void ncio_dense(
    NcIO &ncio,
    VectorMultivec const &vecs,    
    VarSet const &contract,
    ibmisc::Indexing const &indexing,
    std::string const &vname_base = "")
{
    if (ncio.rw != 'w') (*icebin_error)(-1,
        "ncio_dense(VectorMultivec) only writes, no read.");

    // Create the variables
    NcDimSpec dim_spec({"time"}, {-1});
    append(dim_spec, indexing);
    contract.ncdefine(ncio, dim_spec.to_dims(ncio), vname_base);

    ncio += std::bind(ncwrite_dense,
        ncio.nc, &vecs, &contract, &indexing, vname_base);
}
// ------------------------------------------------------------
/** Top-level ncio() to log output from coupler. (coupler->GCM) */
void GCMCoupler::ncio_gcm_input(NcIO &ncio,
    GCMInput &out,        // All MPI ranks included here
    std::array<double,2> &timespan,    // timespan[1] = current time_s
    std::string const &time_units,
    std::string const &vname_base)
{
    ibmisc::ncio_timespan(ncio, timespan, time_units, vname_base);

    for (int iAE=0; iAE<GridAE::count; ++iAE) {
        ncio_dense(ncio, out.gcm_ivalsAE[iAE], gcm_inputsAE[iAE],
            gcm_regridder.indexing(iAE), vname_base);
    }

    ncio_spsparse(ncio, out.E1vE0_s, false, vname_base+"Ev1Ev0");
    ncio_spsparse(ncio, out.AvE1_s, false, vname_base+"AvE1");
    ncio_spsparse(ncio, out.wAvE1_s, false, vname_base+"wAvE1");
    ncio_spsparse(ncio, out.elevE1_s, false, vname_base+"elevE1");
}

/** Top-level ncio() to log input to coupler. (GCM->coupler) */
void GCMCoupler::ncio_gcm_output(NcIO &ncio,
    VectorMultivec const &gcm_ovalsE,
    std::array<double,2> &timespan,    // timespan[1] = current time_s
    std::string const &time_units,
    std::string const &vname_base)
{
    ncio_timespan(ncio, timespan, time_units, vname_base);

    ncio_dense(ncio, gcm_ovalsE, gcm_outputsE,
        gcm_regridder.indexingE, vname_base);
}
// ------------------------------------------------------------
// ======================================================================

}       // namespace
