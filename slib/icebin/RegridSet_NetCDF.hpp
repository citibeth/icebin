#ifndef ICEBIN_REGRID_SET_NETCDF_HPP
#define ICEBIN_REGRID_SET_NETCDF_HPP

namespace icebin {

/** Reads lintransform::Weighted_Compressed out of a NetCDF file */ 
class RegridSet_NetCDF : public RegridSet {

    NcIO ncio;

public:
    std::unique_ptr<ibmisc::lintransform::Weighted> matrix(
        std::string const &spec_name,
        RegridParams const &_params) const;


};

RegridSet_NetCDF::RegridSet_NetCDF(std::string const &fname)
    : ncio(fname, 'r')
{
    // Open the NetCDF file, make sure the Params match
    ...TODO...
}

std::unique_ptr<ibmisc::lintransform::Weighted> RegridSet_NetCDF::matrix(
    std::string const &spec_name) const
{
    std::unique_ptr<lintransform::Weighted_Compressed> M(
        new lintransform::Weighted_Compressed);
    M->ncio(ncio, spec_name);
    return std::unique_ptr<ibmisc::lintransform::Weighted>(M.release());
}



};

#endif
