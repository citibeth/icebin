/*
 * IceBin: A Coupling Library for Ice Models and GCMs
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <mpi.h>        // Must be first

#include <icebin/contracts/contracts.hpp>
#include <icebin/GCMCoupler.hpp>
#include <functional>

namespace icebin {
namespace contracts{

// ======================================================
struct VtableEntry {
    std::function<void(GCMCoupler const &, IceCoupler &)> setup;
};
struct Vtable : public std::map<
    std::pair<GCMCoupler::Type, IceCoupler::Type>,
    VtableEntry >
{
    Vtable();
};


#if defined(BUILD_MODELE) && defined(USE_PISM)
    extern void setup_modele_pism(GCMCoupler const &, IceCoupler &);
#endif

Vtable::Vtable()
{
    VtableEntry entry;

#if defined(BUILD_MODELE) && defined(USE_PISM)
    entry.setup = &setup_modele_pism;
    insert(std::make_pair(
        std::make_pair(GCMCoupler::Type::MODELE, IceCoupler::Type::PISM),
        std::move(entry)));
#endif
}
// -------------------------------------------
static Vtable vtable;

void setup(GCMCoupler const &coupler, IceCoupler &ice_model)
{
    // Dispatch based on BOTH types.
    vtable.at(std::make_pair(coupler.type, ice_model.type))
        .setup(coupler, ice_model);
}
// ======================================================



// -----------------------------------------------------------

std::string flags_to_str(unsigned int flags)
{
    std::string ret = "";
    if (flags & INITIAL) ret += "INITIAL|";

    if (ret.size() > 0) ret.resize(ret.size()-1);
    return ret;
}

}}

