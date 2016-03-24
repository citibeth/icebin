#include <mpi.h>        // Must be first

#include <icebin/contracts/contracts.hpp>
#include <icebin/GCMCoupler.hpp>
#include <functional>

namespace icebin {
namespace contracts{

// ======================================================
struct VtableEntry {
	std::function<void(GCMCoupler &, IceModel &)> setup;
};
struct Vtable : public std::map<
	std::pair<GCMCoupler::Type, IceModel::Type>,
	VtableEntry >
{
	Vtable();
};


#if defined(USE_MODELE) && defined(USE_PISM)
	extern setup_modele_pism(GCMCoupler &, IceModel &);
#endif

Vtable::Vtable()
{
	VtableEntry entry;

#if defined(USE_MODELE) && defined(USE_PISM)
	entry.setup = &setup_modele_pism;
	vtable.insert(std::make_pair(
		std::make_pair(GCMCoupler::Type::MODELE, IceModel::Type::PISM),
		std::move(entry)));
#endif
}
// -------------------------------------------
static Vtable vtable;

void setup(GCMCoupler &coupler, IceModel &ice_model)
{
	vtable.at(std::make_pair(coupler.type, ice_model.type))
		.setup(coupler, ice_model);


	// Check the contract for errors
	VarSet const &ocontract(ice_model.contract[IceModel::OUTPUT]);
    for (size_t i=0; i < ocontract.index.size(); ++i) {
        VarMeta const &cf = ocontract.data[i];

		if ((cf.flags & contracts::GRID_BITS) == contracts::ICE) (*icebin_error)(-1,
			"ERROR: Ice model outputs must be all on the ice grid, field %s is not", cf.name.c_str());
	}


}
// ======================================================



// -----------------------------------------------------------

std::string flags_to_str(unsigned int flags)
{
	std::string ret = "";
	switch(flags & GRID_BITS) {
		case ATMOSPHERE :
			ret += "ATMOSPHERE|";
			break;
		case ICE:
			ret += "ICE|";
			break;
		case ELEVATION:
			ret += "ELEVATION|";
			break;
		default: ;
	}

	if (flags & INITIAL) ret += "INITIAL|";

	if (ret.size() > 0) ret.resize(ret.size()-1);
	return ret;
}

}}

