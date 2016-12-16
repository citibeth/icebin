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
#include <icebin/error.hpp>
#include <ibmisc/indexing.hpp>

namespace icebin {

/** A segment of elevation classes (see add_fhc.py) */
struct HCSegmentData {
    std::string name;
    int base;    // First elevation class of this segment
    int size;    // Number of elevation classes in this segment

    HCSegmentData(std::string const &_name, int _base, int _size)
        : name(_name), base(_base), size(_size) {}
};

/** Parameters passed from the GCM through to the ice model.
These parameters cannot be specific to either the ice model or the GCM.
TODO: Make procedure to read rundeck params and set this stuff up. */
struct GCMParams {
    // ------- Passed into GCMCoupler::allocate()
    ibmisc::Domain domainA, domainA_global;
    MPI_Comm gcm_comm;
    int gcm_root;

    std::string icebin_config_fname;
    std::string config_dir; // Where to look for Ice Model configuration files
    std::string run_dir;    // The GCM run directory

    bool icebin_logging = true ;    // Should IceBin log input & output?

    // Should IceBin update topography?
    bool dynamic_topo = false;

    std::vector<HCSegmentData> hc_segments {
        HCSegmentData("legacy", 0, 1),
        HCSegmentData("sealand", 1, 2),
        HCSegmentData("ec", 3, -1)};    // Last segment must be called ec

    HCSegmentData &ec_segment()
    {
        auto &ec(hc_segments[hc_segments.size()-1]);
        if (ec.name != "ec") (*icebin_error)(-1,
            "The last elevation class segment must be called 'ec'");
        return ec;
    }
};


}
