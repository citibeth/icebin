#ifndef ICEBIN_MODELE_GLOBAL_EC_HPP
#define ICEBIN_MODELE_GLOBAL_EC_HPP

#include <ibmisc/netcdf.hpp>
#include <icebin/modele/hntr.hpp>
#include <icebin/GridSpec.hpp>

namespace icebin {
namespace modele {
namespace global_ec {

/** Metadata read from the file, which we transfer the the output */
struct Metadata {
    bool mismatched;
    icebin::HntrSpec hspecA, hspecI, hspecI2;
    ibmisc::Indexing indexingI, indexingI2, indexingA, indexingHC, indexingE;
    std::vector<double> hcdefs;

    void ncio(ibmisc::NcIO &ncio);
};

inline void Metadata::ncio(ibmisc::NcIO &ncio)
{
    ibmisc::get_or_put_att(*ncio.nc, ncio.rw, "mismatched", mismatched);

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
