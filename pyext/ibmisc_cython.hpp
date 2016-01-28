#pragma once

#include <ibmisc/netcdf.hpp>
#include <ibmisc/ibmisc.hpp>

namespace ibmisc {

static inline ibmisc::NcIO *new_ncio(std::string fname, std::string sfileMode)
{
	netCDF::NcFile::FileMode fileMode;
	if (sfileMode == "read") fileMode = netCDF::NcFile::FileMode::read;
	else if (sfileMode == "write") fileMode = netCDF::NcFile::FileMode::write;
	else if (sfileMode == "replace") fileMode = netCDF::NcFile::FileMode::replace;
	else if (sfileMode == "newFile") fileMode = netCDF::NcFile::FileMode::newFile;
	else {
		(*ibmisc::ibmisc_error)(-1,
			"Bad file mode for netCDF::NcFile::FileMode: %s", sfileMode.c_str());
	}


	return new ibmisc::NcIO(fname, fileMode);
}

}
