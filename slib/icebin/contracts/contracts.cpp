#include <glint2/contracts/contracts.hpp>
#include <functional>

namespace icebin {
namespace contract{

// ======================================================
struct VtableEntry {
	std::function<void(GCMCoupler *, IceModel *)> setup;
};
struct Vtable : public std::map<
	std::pair<GCMCoupler::Type, IceModel::Type>,
	VtableEntry >
{
	Vtable();
};


#if defined(USE_MODELE) && defined(USE_PISM)
	extern setup_modele_pism(GCMCoupler *, IceModel *);
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

void setup(GCMCoupler *coupler, IceModel *model)
{
	vtable.at(coupler.type, model.type)
		.setup(coupler, model);


	// Check the contract for errors
	giss::CouplingContract const &ocontract(ice_model->contract[IceModel::OUTPUT]);
	int nfields = ocontract.size_nounit();
	for (int i=0; i<nfields; ++i) {
		giss::CoupledField const &cf(ocontract.field(i));

		if ((cf.flags & contracts::GRID_BITS) == contracts::ICE) {
			fprintf(stderr, "ERROR: Ice model outputs must be all on the ice grid, field %s is not\n", cf.name.c_str());
			giss::exit(1);
		}
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

