#include <mpi.h>	// Must be first
#include <glint2/IceModel.hpp>
#include <glint2/GCMCoupler.hpp>

namespace glint2 {

IceModel::IceModel(IceModel::Type _type, std::string const &_name, GCMCoupler const *_coupler)
	: type(_type), name(_name), coupler(_coupler), ice_constants(&_coupler->ut_system),
	contract({giss::CouplingContract(), giss::CouplingContract()})
 {}

IceModel::~IceModel() {}


giss::CouplingContract *IceModel::new_CouplingContract() {
	_extra_contracts.push_back(
		std::unique_ptr<giss::CouplingContract>(
		new giss::CouplingContract()));
	return _extra_contracts[_extra_contracts.size()-1].get();
}

}
