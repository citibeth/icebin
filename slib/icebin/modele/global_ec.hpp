#ifndef ICEBIN_MODELE_GLOBAL_EC_HPP
#define ICEBIN_MODELE_GLOBAL_EC_HPP

/** Stuff used to read the output of global_ec.cpp */

#include <ibmisc/netcdf.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/GridSpec.hpp>

namespace icebin {
namespace modele {

BOOST_ENUM_VALUES(GCMGridOption, int,
    (atmosphere) (0)
    (ocean) (1)
    (mismatched) (2)    // (mismatched atmosphere)
)

namespace global_ec {

/** Metadata read from the file, which we transfer the the output */
struct Metadata {
    double eq_rad;    // Radius ouf the earth; must be the same as in ModelE
    GCMGridOption gcm_grid_option;
    icebin::HntrSpec hspecA, hspecI, hspecI2;
    ibmisc::Indexing indexingI, indexingI2, indexingA, indexingHC, indexingE;
    std::vector<double> hcdefs;
    std::vector<int16_t> underice_hc;

    void ncio(ibmisc::NcIO &ncio);
};

inline void Metadata::ncio(ibmisc::NcIO &ncio)
{
    ibmisc::get_or_put_att_enum(*ncio.nc, ncio.rw, "gcm_grid_option", gcm_grid_option);
    ibmisc::get_or_put_att(*ncio.nc, ncio.rw, "eq_rad", "double", &eq_rad, 1);

    std::string const &hspecA_name
        (gcm_grid_option == GCMGridOption::ocean ? "hspecO" : "hspecA");
    std::string const &indexingA_name
        (gcm_grid_option == GCMGridOption::ocean ? "indexingO" : "indexingA");

    hspecA.ncio(ncio, hspecA_name);
    hspecI.ncio(ncio, "hspecI");
    hspecI2.ncio(ncio, "hspecI2");

    indexingI.ncio(ncio, "indexingI");
    indexingI2.ncio(ncio, "indexingI2");
    indexingA.ncio(ncio, indexingA_name);
    indexingHC.ncio(ncio, "indexingHC");
    indexingE.ncio(ncio, "indexingE");

    auto _nhc(get_or_add_dims(ncio, {"nhc"}, {hcdefs.size()}));
    ncio_vector(ncio, hcdefs, true, "hcdefs", "double", _nhc);
    ncio_vector(ncio, underice_hc, true, "underice_hc", "short", _nhc);  // Must be short for NetCDF3
}

}}}    // namespace
#endif    // guard
