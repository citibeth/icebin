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

#pragma once

#include <boost/filesystem.hpp>
#include <ibmisc/time.hpp>

namespace icebin {

/** Parameters passed from the GCM through to the ice model.
These parameters cannot be specific to either the ice model or the GCM. */
struct GCMParams {
    boost::filesystem::path config_dir; // Where to look for Ice Model configuration files
    boost::filesystem::path run_dir;
    ibmisc::time::tm time_base; // Corresponds to time_s == 0
    std::string time_units;     // CF-compliant string describing the time units
    double time_start_s;        // Start of simulation, as far as ice model is concerned (seconds since time_base).

    GCMParams();

    GCMParams(
        boost::filesystem::path const &_config_dir);

    void set_start_time(
        ibmisc::time::tm const &_time_base,
        double time_start_s);
};


}
