#ifndef ICEBIN_ELEVMASK_HPP
#define ICEBIN_ELEVMASK_HPP

#include <memory>
#include <blitz/array.h>
#include <ibmisc/netcdf.hpp>

namespace icebin {


/** Reads and allocate
@param emI Elevation-mask for continental area (ice and bare land)
@param emI_ice Elevation-mask for just ice-covered areas */
void read_elevmask_pism(
    std::string const &fname,
    int const itime,    // Value of PISM time dimension to read
    blitz::Array<double,1> &emI_land,
    blitz::Array<double,1> &emI_ice);


}    // namespace
#endif    // guard
