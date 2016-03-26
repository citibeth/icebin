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
#include <boost/format.hpp>
#include <icebin/GCMParams.hpp>

namespace icebin {

GCMParams::GCMParams() :
    gcm_rank(-1), gcm_root(-1)
{
    // Set a reasonable default initial value for start time.
    ibmisc::time::tm const time_base(1900,1,1);
    set_start_time(time_base, -1);
}

void GCMParams::set_start_time(
    ibmisc::time::tm const &_time_base,
    double _time_start_s)
{
    printf("GCMParams::set_start_time: time_units = %s\n", time_units.c_str());
    time_base = _time_base;
    time_start_s = _time_start_s;
    time_units = str(boost::format("seconds since %04d-%02d-%02d %02d:%02d:%02d")
        % time_base.year()
        % time_base.month()
        % time_base.mday()
        % time_base.tm_hour
        % time_base.tm_min
        % time_base.tm_sec);
    printf("GCMParams::set_start_time: %s\n", time_units.c_str());
}

GCMParams::GCMParams(
    MPI_Comm const _gcm_comm,
    int _gcm_root,
    boost::filesystem::path const &_config_dir)
: gcm_comm(_gcm_comm), gcm_root(_gcm_root), config_dir(_config_dir), time_start_s(0)
{
    MPI_Comm_rank(gcm_comm, &gcm_rank);

}

}

