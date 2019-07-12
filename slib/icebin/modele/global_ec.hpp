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

    void ncio(ibmisc::NcIO &ncio);
};

inline void Metadata::ncio(ibmisc::NcIO &ncio)
{
    ibmisc::get_or_put_att_enum(*ncio.nc, ncio.rw, "gcm_grid_option", gcm_grid_option);
    ibmisc::get_or_put_att(*ncio.nc, ncio.rw, "eq_rad", eq_rad);

    hspecA.ncio(ncio, "hspecA");
    hspecI.ncio(ncio, "hspecI");
    hspecI2.ncio(ncio, "hspecI2");

    indexingI.ncio(ncio, "indexingI");
    indexingI2.ncio(ncio, "indexingI2");
    indexingA.ncio(ncio, "indexingA");
    indexingHC.ncio(ncio, "indexingHC");
    indexingE.ncio(ncio, "indexingE");

    ncio_vector(ncio, hcdefs, true, "hcdefs", "double",
        get_or_add_dims(ncio, {"nhc"}, {hcdefs.size()}));

}

}}}    // namespace
#endif    // guard
