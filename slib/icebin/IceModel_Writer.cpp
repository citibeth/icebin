/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <mpi.h>        // Must be first
#include <boost/filesystem.hpp>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <icebin/GCMParams.hpp>
#include <icebin/GCMCoupler.hpp>
#include <ibmisc/string.hpp>
#include <ibmisc/netcdf.hpp>

using namespace ibmisc;
using namespace netCDF;


namespace icebin {

/** Specialized init signature for IceModel_Writer */
void IceModel_Writer::init(
    IceModel::IO _io,   // Ic
    IceModel const *_main_model)

{
    printf("BEGIN IceModel_Writer::init(%s)\n", name().c_str());

    main_model = _main_model;
    io = _io;
    output_file_initialized = false;

    // Try to be clever about making multi-dimensional arrays
    // in the output according to the grid the user expects.
    auto &indexing(gridI()->indexing);
    dim_names = {"time"};
    counts = {1};
    cur = {0};
    for (size_t i=0; i<gridI()->indexing.rank(); ++i) {
        // ix goes 0...n-1 for row-major, n-1..0 for col-major
        int ix = indexing.indices[i];
        dim_names.push_back(string_printf("dim%d", ix));
        counts.push_back(indexing.extent[ix]);
        cur.push_back(0);
    }

    // Only need to run one copy of this
    GCMParams const &gcm_params(coupler->gcm_params);
    if (gcm_params.gcm_rank != gcm_params.gcm_root) return;

    // Put our output files in this directory, one named per ice sheet.
    auto output_dir = boost::filesystem::absolute(
        boost::filesystem::path(
            (io == IceModel::INPUT ? "ice_model_in" : "ice_model_out")),
        coupler->gcm_params.run_dir);
    boost::filesystem::create_directory(output_dir);    // Make sure it exists

    // Set up the output file
    // Create netCDF variables based on details of the coupling contract.xs
    output_fname = (output_dir / (name() + ".nc")).string();
    printf("IceModel_Writer opening file %s\n", output_fname.c_str());

    printf("END IceModel_Writer::init_from_ice_model(%s)\n", name().c_str());
}

void IceModel_Writer::start_time_set()
{
    // We just need the input (or output) coupling contract
    contract[io] = main_model->contract[io];
}

void IceModel_Writer::init_output_file()
{
    GCMParams const &gcm_params(coupler->gcm_params);

    printf("BEGIN IceModel_Writer::init_output_file(%s)\n", output_fname.c_str());

    NcIO ncio(output_fname, NcFile::replace);

    auto dims(get_or_add_dims(ncio, dim_names, counts));

    NcDim one_dim = ncio.nc->addDim("one", 1);
    NcVar info_var = ncio.nc->addVar("grid", ibmisc::get_nc_type<double>(), one_dim);

    info_var.putAtt("file", coupler->fname);
    info_var.putAtt("variable", coupler->vname);
    info_var.putAtt("ice_sheet", this->name());

    NcVar time0_var = ncio.nc->addVar("time0", ibmisc::get_nc_type<double>(), one_dim);
    time0_var.putAtt("units", gcm_params.time_units);
    time0_var.putAtt("calendar", "365_day");
    time0_var.putAtt("axis", "T");
    time0_var.putAtt("long_name", "Simulation start time");

    NcVar time_var = ncio.nc->addVar("time", ibmisc::get_nc_type<double>(), dims[0]);
    time_var.putAtt("units", gcm_params.time_units);
    time_var.putAtt("calendar", "365_day");
    time_var.putAtt("axis", "T");
    time_var.putAtt("long_name", "Coupling times");


    for (size_t i=0; i < contract[io].index.size(); ++i) {
        VarMeta &cf = contract[io].data[i];
        NcVar var = ncio.nc->addVar(cf.name, ibmisc::get_nc_type<double>(), dims);
        var.putAtt("units", cf.units);
        var.putAtt("description", cf.description);
    }

    // Put initial time in it...
    time0_var.putVar({0}, {1}, &gcm_params.time_start_s);
#if 0
    std::vector<size_t> cur_b {0};
    std::vector<size_t> counts_b {1};
    long cur_b[1]{0};
    long counts_b[1]{1};
    time0_var.set_cur(cur_b);
    time0_var.put(&gcm_params.time_start_s, counts_b);
#endif
    ncio.close();

    output_file_initialized = true;

    printf("END IceModelWriter::init_output_file(%s)\n", output_fname.c_str());
}

/** @param index Index of each grid value.
@param vals The values themselves -- could be SMB, Energy, something else...
TODO: More params need to be added.  Time, return values, etc. */
void IceModel_Writer::run_decoded(double time_s,
    std::vector<blitz::Array<double,1>> const &ivalsI)
{
    // Only need to run one copy of this
    GCMParams const &gcm_params(coupler->gcm_params);
    if (gcm_params.gcm_rank != gcm_params.gcm_root) return;

printf("BEGIN IceModel_Writer::run_decoded\n");
    if (!output_file_initialized) init_output_file();
    NcFile nc(output_fname.c_str(), NcFile::write);

    // Read index info
    cur[0] = nc.getDim("time").getSize();

    // Write the current time
    NcVar time_var = nc.getVar("time");
    time_var.putVar(cur, counts, &time_s);


    // Write the other variables
    for (size_t i=0; i < contract[io].index.size(); ++i) {
        VarMeta &cf = contract[io].data[i];

        if (cf.name == "unit") continue;

        NcVar ncvar = nc.getVar(cf.name.c_str());
        ncvar.putVar(cur, counts, ivalsI[i].data());
    }

    nc.close();
printf("END IceModel_Writer::run_decoded\n");
    
}


}
