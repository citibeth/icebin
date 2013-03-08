#include <glint2/GCMCoupler.hpp>

namespace glint2 {

/** Query all the ice models to figure out what fields they need */
std::set<IceField> GCMCoupler::get_required_fields()
{
	std::set<IceField> ret;
	for (auto model = models.begin(); model != models.end(); ++model) {
		model->get_required_fields(ret);
	}
	return ret;
}


void GCMCoupler::read_from_netcdf(NcFile &nc, std::string const &vname,
	std::vector<string> const &sheet_names)
{
	int i = 0;
	for (auto name = sheet_names.begin(); name != sheet_names.end(); ++name) {
		models.insert(i, read_icemodel(nc, vname + "." + *name));
		++i;
	}
}


}