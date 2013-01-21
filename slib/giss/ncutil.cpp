#include <giss/ncutil.hpp>

namespace giss {

NcDim *get_or_add_dim(NcFile &nc, std::string const &dim_name, long dim_size)
{
	// Look up dim the slow way...
	int num_dims = nc.num_dims();
	for (int i=0; i<num_dims; ++i) {
		NcDim *dim = nc.get_dim(i);
		char const *name = dim->name();

		if (strcmp(name, dim_name.c_str()) == 0) {
			long sz = dim->size();
			if (sz != dim_size) {
				fprintf(stderr, "Error: dimension %s (size = %ld) being redefined to size = %ld\n",
					dim_name.c_str(), sz, dim_size);
				throw std::exception();
			}
			return dim;
		}
	}

	return nc.add_dim(dim_name.c_str(), dim_size);
}

NcVar *get_var_safe(NcFile &nc, std::string const &var_name)
{
	// Look up var the slow way...
	int num_vars = nc.num_vars();
	for (int i=0; i<num_vars; ++i) {
		NcVar *var = nc.get_var(i);
		char const *name = var->name();

		if (strcmp(name, var_name.c_str()) == 0) {
			return var;
		}
	}
	return NULL;
}

std::vector<double> read_double_vector(NcFile &nc, std::string const &var_name)
{
	// Read points vector
	NcVar *vpoints = nc.get_var(var_name.c_str());
	long npoints = vpoints->get_dim(0)->size();
	std::vector<double> points(npoints);
	vpoints->get(&points[0], npoints);
	return points;
}

std::vector<int> read_int_vector(NcFile &nc, std::string const &var_name)
{
	// Read points vector
	NcVar *vpoints = nc.get_var(var_name.c_str());
	long npoints = vpoints->get_dim(0)->size();
	std::vector<int> points(npoints);
	vpoints->get(&points[0], npoints);
	return points;
}


}
