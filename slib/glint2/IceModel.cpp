#include <mpi.h>	// Must be first
#include <glint2/IceModel.hpp>
#include <glint2/GCMCoupler.hpp>

namespace glint2 {

IceModel::IceModel(IceModel::Type _type, GCMCoupler const *_coupler)
	: type(_type), coupler(_coupler), ice_constants(&_coupler->ut_system),
	contract({giss::CouplingContract(&coupler->ut_system), giss::CouplingContract(&coupler->ut_system)})
 {}

IceModel::~IceModel() {}


giss::CouplingContract *IceModel::new_CouplingContract() {
	_extra_contracts.push_back(
		std::unique_ptr<giss::CouplingContract>(
		new giss::CouplingContract(&coupler->ut_system)));
	return _extra_contracts[_extra_contracts.size()-1].get();
}

}
