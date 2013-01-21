#pragma once

#include <cstring>
#include <vector>
#include <netcdfcpp.h>

namespace giss {

extern NcDim *get_or_add_dim(NcFile &nc, std::string const &dim_name, long dim_size);

NcVar *get_var_safe(NcFile &nc, std::string const &var_name);

extern std::vector<double> read_double_vector(NcFile &nc, std::string const &var_name);

extern std::vector<int> read_int_vector(NcFile &nc, std::string const &var_name);


}
